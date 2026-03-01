// ═══════════════════════════════════════════════════════════════════════════════
// EMBER CherchiBackend — Single-Pass Shatter + Open-Mesh Plane Cutter (v2.0)
// ═══════════════════════════════════════════════════════════════════════════════
//
// FIXES IN v2.0:
//   1. Open mesh B (plane/grid) detected via boundary edges → plane-side
//      classification instead of winding number (which is always 0 for open).
//   2. Minimum-area shard culling after split (1e-12 threshold).
//   3. Houdini Shatter convention: B is cutter → 2 pieces, not 3.

#include "CherchiBackend.h"
#include "../ember/PolygonSoup.h"
#include "../ember/ExactPredicates.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>

namespace ember {

// ════════════════════════════════════════════════════════════════════════════
// Constants
// ════════════════════════════════════════════════════════════════════════════

static constexpr uint32_t CAP_FLAG       = 0x80000000u;
static constexpr double   MIN_TRI_AREA   = 1e-12;  // shard culling threshold

// ════════════════════════════════════════════════════════════════════════════
// Forward declarations
// ════════════════════════════════════════════════════════════════════════════

struct CherchiBackend::ImplicitPoint    { uint32_t type; uint32_t data[4]; };
struct CherchiBackend::ArrangementPoint { uint32_t index; bool is_explicit; };
struct CherchiBackend::TriangleIntersection { uint32_t triA, triB, num_points; };
struct CherchiBackend::LocalArrangement { std::vector<uint32_t> triangles, vertices; };

// ════════════════════════════════════════════════════════════════════════════
// §1  Double-precision geometry
// ════════════════════════════════════════════════════════════════════════════

using D3 = std::array<double, 3>;

struct DTri {
    D3       v[3];
    int      mesh_id;
    uint32_t src_prim;
};

static D3 d3sub(const D3& a, const D3& b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}
static D3 d3cross(const D3& a, const D3& b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
static double d3dot(const D3& a, const D3& b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
static double d3len(const D3& a) {
    return std::sqrt(d3dot(a,a));
}
static D3 d3_lerp(const D3& a, const D3& b, double t) {
    return {a[0]+t*(b[0]-a[0]), a[1]+t*(b[1]-a[1]), a[2]+t*(b[2]-a[2])};
}
static D3 d3_centroid(const DTri& t) {
    return {(t.v[0][0]+t.v[1][0]+t.v[2][0])/3.0,
            (t.v[0][1]+t.v[1][1]+t.v[2][1])/3.0,
            (t.v[0][2]+t.v[1][2]+t.v[2][2])/3.0};
}
static D3 triToD3(const Triangle& t, int vi) {
    return {(double)t.v[vi][0], (double)t.v[vi][1], (double)t.v[vi][2]};
}
static DTri toDTri(const Triangle& t) {
    DTri d;
    for (int i=0;i<3;++i) d.v[i]=triToD3(t,i);
    d.mesh_id=t.mesh_id; d.src_prim=t.src_prim_idx;
    return d;
}
static Triangle fromDTri(const DTri& d, int piece_id, bool is_cap) {
    Triangle t;
    for (int i=0;i<3;++i) {
        t.v[i][0]=(float)d.v[i][0]; t.v[i][1]=(float)d.v[i][1]; t.v[i][2]=(float)d.v[i][2];
        t.iv[i]={0,0,0};
    }
    t.mesh_id = piece_id;
    t.src_prim_idx = is_cap ? (d.src_prim | CAP_FLAG) : d.src_prim;
    return t;
}
static DTri flipDTri(DTri d) { std::swap(d.v[1], d.v[2]); return d; }

// Triangle area (2x cross product magnitude)
static double triArea(const DTri& t) {
    D3 e0 = d3sub(t.v[1], t.v[0]);
    D3 e1 = d3sub(t.v[2], t.v[0]);
    D3 cr = d3cross(e0, e1);
    return 0.5 * d3len(cr);
}

// ════════════════════════════════════════════════════════════════════════════
// §2  Exact predicates
// ════════════════════════════════════════════════════════════════════════════

static int exactSign(const D3& A, const D3& B, const D3& C, const D3& P) {
    return predicates::orient3d_filtered(
        A[0],A[1],A[2], B[0],B[1],B[2], C[0],C[1],C[2], P[0],P[1],P[2]);
}
static double planeDistD(const D3& A, const D3& B, const D3& C, const D3& P) {
    D3 ab=d3sub(B,A), ac=d3sub(C,A);
    D3 n=d3cross(ab,ac);
    return d3dot(n,d3sub(P,A));
}

// ════════════════════════════════════════════════════════════════════════════
// §3  Open mesh detection (boundary edge test)
// ════════════════════════════════════════════════════════════════════════════

struct HalfEdge {
    uint32_t a, b;  // Uses point indices from a local index map
};
struct HalfEdgeHash {
    size_t operator()(const std::pair<uint32_t,uint32_t>& e) const {
        return std::hash<uint64_t>()((uint64_t)e.first << 32 | e.second);
    }
};

// Returns true if mesh B has boundary edges (open surface / plane / grid)
static bool isMeshOpen(const std::vector<DTri>& mesh) {
    // Count directed half-edges. A manifold closed mesh has every
    // half-edge paired with its reverse. An open mesh has unpaired edges.
    // We use snapped vertex indices for matching.
    
    std::unordered_map<uint64_t, int> edgeCount;
    auto snapKey = [](const D3& v) -> uint64_t {
        int32_t x = (int32_t)std::round(v[0]*1e6);
        int32_t y = (int32_t)std::round(v[1]*1e6);
        int32_t z = (int32_t)std::round(v[2]*1e6);
        return ((uint64_t)(uint32_t)x << 40) ^ ((uint64_t)(uint32_t)y << 20) ^ (uint32_t)z;
    };
    
    for (const auto& tri : mesh) {
        uint64_t vi[3];
        for (int i=0;i<3;++i) vi[i] = snapKey(tri.v[i]);
        for (int e=0;e<3;++e) {
            uint64_t a = vi[e], b = vi[(e+1)%3];
            uint64_t fwd = (a < b) ? (a*2654435761ULL + b) : (b*2654435761ULL + a);
            edgeCount[fwd]++;
        }
    }
    
    for (const auto& kv : edgeCount) {
        if (kv.second != 2) return true;  // boundary edge (shared by 1 or >2 tris)
    }
    return false;
}

// ════════════════════════════════════════════════════════════════════════════
// §4  Best-fit plane for open mesh B
// ════════════════════════════════════════════════════════════════════════════

struct FitPlane {
    D3 origin;
    D3 normal;
    bool valid;
};

static FitPlane fitPlaneToMesh(const std::vector<DTri>& mesh) {
    FitPlane fp;
    fp.valid = false;
    if (mesh.empty()) return fp;
    
    // Area-weighted average normal + centroid
    D3 centroid = {0,0,0};
    D3 normal   = {0,0,0};
    double totalArea = 0;
    int count = 0;
    
    for (const auto& tri : mesh) {
        D3 e0 = d3sub(tri.v[1], tri.v[0]);
        D3 e1 = d3sub(tri.v[2], tri.v[0]);
        D3 n  = d3cross(e0, e1);
        double area = 0.5 * d3len(n);
        
        normal[0] += n[0]; normal[1] += n[1]; normal[2] += n[2];
        for (int v=0;v<3;++v) {
            centroid[0] += tri.v[v][0];
            centroid[1] += tri.v[v][1];
            centroid[2] += tri.v[v][2];
        }
        totalArea += area;
        count += 3;
    }
    
    if (count == 0 || totalArea < 1e-30) return fp;
    
    centroid[0] /= count; centroid[1] /= count; centroid[2] /= count;
    double nlen = d3len(normal);
    if (nlen < 1e-30) return fp;
    normal[0] /= nlen; normal[1] /= nlen; normal[2] /= nlen;
    
    fp.origin = centroid;
    fp.normal = normal;
    fp.valid  = true;
    return fp;
}

// Classify point against fitted plane using exact orient3d
// Returns +1 (positive side), -1 (negative side), 0 (on plane)
static int classifyAgainstFitPlane(const D3& P, const FitPlane& fp) {
    // Build three points defining the plane
    D3 A = fp.origin;
    // Construct two tangent vectors orthogonal to normal
    D3 up = {0,1,0};
    if (std::abs(d3dot(fp.normal, up)) > 0.9) up = {1,0,0};
    D3 t1 = d3cross(fp.normal, up);
    double t1len = d3len(t1);
    if (t1len < 1e-30) return 0;
    t1[0]/=t1len; t1[1]/=t1len; t1[2]/=t1len;
    D3 t2 = d3cross(fp.normal, t1);
    
    D3 B = {A[0]+t1[0], A[1]+t1[1], A[2]+t1[2]};
    D3 C = {A[0]+t2[0], A[1]+t2[1], A[2]+t2[2]};
    
    return exactSign(A, B, C, P);
}

// ════════════════════════════════════════════════════════════════════════════
// §5  AABB + tri-tri test
// ════════════════════════════════════════════════════════════════════════════

struct BBox {
    double mn[3]={1e30,1e30,1e30}, mx[3]={-1e30,-1e30,-1e30};
    void expand(const D3& v){for(int i=0;i<3;++i){mn[i]=std::min(mn[i],v[i]);mx[i]=std::max(mx[i],v[i]);}}
    bool overlaps(const BBox& o)const{for(int i=0;i<3;++i)if(mn[i]>o.mx[i]||mx[i]<o.mn[i])return false;return true;}
    static BBox ofTri(const Triangle& t){BBox b;for(int v=0;v<3;++v)b.expand(triToD3(t,v));return b;}
};

static bool triTriIntersect(const Triangle& A, const Triangle& B) {
    D3 a[3],b[3];
    for(int i=0;i<3;++i){a[i]=triToD3(A,i);b[i]=triToD3(B,i);}
    int sA[3],sB[3];
    for(int i=0;i<3;++i)sA[i]=exactSign(b[0],b[1],b[2],a[i]);
    if((sA[0]>0&&sA[1]>0&&sA[2]>0)||(sA[0]<0&&sA[1]<0&&sA[2]<0))return false;
    if(sA[0]==0&&sA[1]==0&&sA[2]==0)return false;
    for(int i=0;i<3;++i)sB[i]=exactSign(a[0],a[1],a[2],b[i]);
    if((sB[0]>0&&sB[1]>0&&sB[2]>0)||(sB[0]<0&&sB[1]<0&&sB[2]<0))return false;
    if(sB[0]==0&&sB[1]==0&&sB[2]==0)return false;
    return true;
}

// ════════════════════════════════════════════════════════════════════════════
// §6  Exact-sign plane splitting
// ════════════════════════════════════════════════════════════════════════════

static void splitByPlaneExact(
    const DTri& T, const D3& pA, const D3& pB, const D3& pC,
    std::vector<DTri>& pos, std::vector<DTri>& neg)
{
    int s[3]; double d[3];
    for(int i=0;i<3;++i){
        s[i]=exactSign(pA,pB,pC,T.v[i]);
        d[i]=planeDistD(pA,pB,pC,T.v[i]);
    }
    int npos=(s[0]>0)+(s[1]>0)+(s[2]>0);
    int nneg=(s[0]<0)+(s[1]<0)+(s[2]<0);
    int nzer=(s[0]==0)+(s[1]==0)+(s[2]==0);

    if(nneg==0){pos.push_back(T);return;}
    if(npos==0){neg.push_back(T);return;}

    if(nzer==1){
        int z=-1; for(int i=0;i<3;++i)if(s[i]==0)z=i;
        int i1=(z+1)%3, i2=(z+2)%3;
        if(s[i1]==s[i2]){if(s[i1]>0)pos.push_back(T);else neg.push_back(T);return;}
        double denom=d[i1]-d[i2];
        double t=(std::abs(denom)>1e-30)?(d[i1]/denom):0.5;
        t=std::max(0.0,std::min(1.0,t));
        D3 P=d3_lerp(T.v[i1],T.v[i2],t);
        DTri t1;t1.mesh_id=T.mesh_id;t1.src_prim=T.src_prim;
        t1.v[0]=T.v[z];t1.v[1]=T.v[i1];t1.v[2]=P;
        DTri t2;t2.mesh_id=T.mesh_id;t2.src_prim=T.src_prim;
        t2.v[0]=T.v[z];t2.v[1]=P;t2.v[2]=T.v[i2];
        if(s[i1]>0){pos.push_back(t1);neg.push_back(t2);}
        else       {neg.push_back(t1);pos.push_back(t2);}
        return;
    }
    if(nzer>=2){if(npos>0)pos.push_back(T);else neg.push_back(T);return;}

    int lone=-1,o0=-1,o1=-1;
    int lone_sign=(npos==1)?1:-1;
    for(int i=0;i<3;++i){if(s[i]==lone_sign)lone=i;else{if(o0<0)o0=i;else o1=i;}}
    if(lone<0||o0<0||o1<0){pos.push_back(T);return;}

    auto safeLerp=[](double dA,double dB)->double{
        double den=dA-dB;if(std::abs(den)<1e-30)return 0.5;
        double t=dA/den;return std::max(0.0,std::min(1.0,t));};
    D3 P0=d3_lerp(T.v[lone],T.v[o0],safeLerp(d[lone],d[o0]));
    D3 P1=d3_lerp(T.v[lone],T.v[o1],safeLerp(d[lone],d[o1]));

    DTri lT;lT.mesh_id=T.mesh_id;lT.src_prim=T.src_prim;
    lT.v[0]=T.v[lone];lT.v[1]=P0;lT.v[2]=P1;
    DTri oA;oA.mesh_id=T.mesh_id;oA.src_prim=T.src_prim;
    oA.v[0]=P0;oA.v[1]=T.v[o0];oA.v[2]=T.v[o1];
    DTri oB;oB.mesh_id=T.mesh_id;oB.src_prim=T.src_prim;
    oB.v[0]=P0;oB.v[1]=T.v[o1];oB.v[2]=P1;

    if(lone_sign>0){pos.push_back(lT);neg.push_back(oA);neg.push_back(oB);}
    else           {neg.push_back(lT);pos.push_back(oA);pos.push_back(oB);}
}

// ════════════════════════════════════════════════════════════════════════════
// §7  Winding number via +X ray (for closed mesh B only)
// ════════════════════════════════════════════════════════════════════════════

static int windingNumber(const D3& O, const std::vector<DTri>& mesh) {
    int w=0;
    for(const auto& T:mesh){
        double oy=O[1],oz=O[2];
        int d1=predicates::orient2d_filtered(T.v[0][1],T.v[0][2],T.v[1][1],T.v[1][2],oy,oz);
        int d2=predicates::orient2d_filtered(T.v[1][1],T.v[1][2],T.v[2][1],T.v[2][2],oy,oz);
        int d3=predicates::orient2d_filtered(T.v[2][1],T.v[2][2],T.v[0][1],T.v[0][2],oy,oz);
        if(!(d1>0&&d2>0&&d3>0)&&!(d1<0&&d2<0&&d3<0))continue;
        double abx=T.v[1][0]-T.v[0][0],aby=T.v[1][1]-T.v[0][1],abz=T.v[1][2]-T.v[0][2];
        double acx=T.v[2][0]-T.v[0][0],acy=T.v[2][1]-T.v[0][1],acz=T.v[2][2]-T.v[0][2];
        double nx=aby*acz-abz*acy;
        if(std::abs(nx)<1e-15)continue;
        double ny=abz*acx-abx*acz, nz=abx*acy-aby*acx;
        double t=-(nx*(O[0]-T.v[0][0])+ny*(O[1]-T.v[0][1])+nz*(O[2]-T.v[0][2]))/nx;
        if(t<=0.0)continue;
        w+=(nx>0.0)?+1:-1;
    }
    return w;
}

// ════════════════════════════════════════════════════════════════════════════
// §8  Standard Boolean inclusion
// ════════════════════════════════════════════════════════════════════════════

static bool shouldInclude(int winding, BooleanOp op, bool from_A, bool& flip) {
    flip=false;
    bool inside=(winding!=0);
    switch(op){
        case BooleanOp::Union:        return !inside;
        case BooleanOp::Intersection: return inside;
        case BooleanOp::DiffAB:
            if(from_A)return !inside; flip=true; return inside;
        case BooleanOp::DiffBA:
            if(!from_A)return !inside; flip=true; return inside;
        default: return !inside;
    }
}

// ════════════════════════════════════════════════════════════════════════════
// §9  Split engine — ONE best-fit cut per A triangle
// ════════════════════════════════════════════════════════════════════════════
//
// KEY FIX: Instead of splitting each A triangle by ALL intersecting B
// planes (creating 2^N fragments with staircase debris), we split each A
// triangle by only the SINGLE B triangle whose plane produces the most
// significant cut. This gives at most 3 fragments per A triangle.

struct SplitResult {
    std::vector<DTri> splitA;
    std::vector<DTri> splitB;
};

// Score how "significant" a plane cut is for a triangle:
// Returns the minimum absolute signed distance of any vertex to the plane.
// Higher = more meaningful cut (cuts through the middle, not a corner).
static double cutSignificance(const DTri& T, const DTri& cutter) {
    double dmin = 1e30;
    for (int i = 0; i < 3; ++i) {
        double d = std::abs(planeDistD(cutter.v[0], cutter.v[1], cutter.v[2], T.v[i]));
        dmin = std::min(dmin, d);
    }
    return dmin;
}

static SplitResult splitBothMeshes(
    const PolygonSoup& soup, uint32_t nA, uint32_t nB,
    const std::vector<DTri>& origA, const std::vector<DTri>& origB)
{
    // Compute mesh scale for relative thresholds
    double totalArea = 0;
    for (uint32_t i = 0; i < nA; ++i) totalArea += triArea(origA[i]);
    double avgArea = (nA > 0) ? totalArea / nA : 1e-6;
    double minArea = avgArea * 1e-4;  // Scale-relative shard threshold

    // Detect intersecting pairs
    std::vector<std::vector<uint32_t>> cutsOnA(nA), cutsOnB(nB);
    {
        std::vector<BBox> bboxB(nB);
        for(uint32_t j=0;j<nB;++j) bboxB[j]=BBox::ofTri(soup.triangles[nA+j]);
        for(uint32_t i=0;i<nA;++i){
            BBox bA=BBox::ofTri(soup.triangles[i]);
            for(uint32_t j=0;j<nB;++j){
                if(!bA.overlaps(bboxB[j]))continue;
                if(triTriIntersect(soup.triangles[i],soup.triangles[nA+j])){
                    cutsOnA[i].push_back(j); cutsOnB[j].push_back(i);
                }
            }
        }
    }

    SplitResult sr;

    // Split mesh A: each triangle by its SINGLE best cutting plane
    sr.splitA.reserve(nA * 2);
    for (uint32_t i = 0; i < nA; ++i) {
        if (cutsOnA[i].empty()) { sr.splitA.push_back(origA[i]); continue; }

        // Find the best cutter: highest significance score
        uint32_t bestJ = cutsOnA[i][0];
        double bestScore = cutSignificance(origA[i], origB[bestJ]);
        for (size_t k = 1; k < cutsOnA[i].size(); ++k) {
            uint32_t j = cutsOnA[i][k];
            double score = cutSignificance(origA[i], origB[j]);
            if (score > bestScore) { bestScore = score; bestJ = j; }
        }

        // Split by the single best plane
        std::vector<DTri> p, n;
        splitByPlaneExact(origA[i], origB[bestJ].v[0], origB[bestJ].v[1], origB[bestJ].v[2], p, n);

        for (const auto& f : p) { if (triArea(f) >= minArea) sr.splitA.push_back(f); }
        for (const auto& f : n) { if (triArea(f) >= minArea) sr.splitA.push_back(f); }
    }

    // Split mesh B: same approach, single best cutter from A
    sr.splitB.reserve(nB * 2);
    for (uint32_t j = 0; j < nB; ++j) {
        if (cutsOnB[j].empty()) { sr.splitB.push_back(origB[j]); continue; }

        uint32_t bestI = cutsOnB[j][0];
        double bestScore = cutSignificance(origB[j], origA[bestI]);
        for (size_t k = 1; k < cutsOnB[j].size(); ++k) {
            uint32_t i = cutsOnB[j][k];
            double score = cutSignificance(origB[j], origA[i]);
            if (score > bestScore) { bestScore = score; bestI = i; }
        }

        std::vector<DTri> p, n;
        splitByPlaneExact(origB[j], origA[bestI].v[0], origA[bestI].v[1], origA[bestI].v[2], p, n);

        for (const auto& f : p) { if (triArea(f) >= minArea) sr.splitB.push_back(f); }
        for (const auto& f : n) { if (triArea(f) >= minArea) sr.splitB.push_back(f); }
    }

    return sr;
}

// ════════════════════════════════════════════════════════════════════════════
// §10  Shatter classification (closed B: winding, open B: plane-side)
// ════════════════════════════════════════════════════════════════════════════

static std::vector<Triangle> classifyShatter(
    const SplitResult& sr,
    const std::vector<DTri>& origA, const std::vector<DTri>& origB,
    bool bIsOpen, const FitPlane& bPlane)
{
    std::vector<Triangle> output;
    output.reserve(sr.splitA.size() + sr.splitB.size());

    // Classify A sub-triangles
    for (const auto& frag : sr.splitA) {
        D3 c = d3_centroid(frag);
        bool inside_B;
        if (bIsOpen) {
            // Open B: classify by which side of B's fitted plane
            int side = classifyAgainstFitPlane(c, bPlane);
            inside_B = (side < 0);  // negative side = "inside" the cut
        } else {
            inside_B = (windingNumber(c, origB) != 0);
        }

        if (!inside_B) {
            output.push_back(fromDTri(frag, 0, false));  // Piece 0
        } else {
            output.push_back(fromDTri(frag, 1, false));  // Piece 1
        }
    }

    // B sub-triangles: only those inside A become caps/shells
    for (const auto& frag : sr.splitB) {
        D3 c = d3_centroid(frag);
        bool inside_A;
        if (bIsOpen) {
            // For open B cutter: B's triangles that overlap A's volume
            // Use A's winding (A is always closed for Shatter to make sense)
            inside_A = (windingNumber(c, origA) != 0);
        } else {
            inside_A = (windingNumber(c, origA) != 0);
        }

        if (inside_A) {
            output.push_back(fromDTri(flipDTri(frag), 0, true));  // Cap of Piece 0
            output.push_back(fromDTri(frag, 1, true));            // Cap of Piece 1
        }
        // B outside A → discarded
    }

    return output;
}

// ════════════════════════════════════════════════════════════════════════════
// §11  Standard Boolean classification (closed B: winding, open B: plane)
// ════════════════════════════════════════════════════════════════════════════

static std::vector<Triangle> classifyStandard(
    const SplitResult& sr,
    const std::vector<DTri>& origA, const std::vector<DTri>& origB,
    BooleanOp op,
    bool bIsOpen, const FitPlane& bPlane)
{
    std::vector<Triangle> output;
    output.reserve(sr.splitA.size() + sr.splitB.size());

    for (const auto& frag : sr.splitA) {
        D3 c = d3_centroid(frag);
        int w;
        if (bIsOpen) {
            int side = classifyAgainstFitPlane(c, bPlane);
            w = (side < 0) ? 1 : 0;
        } else {
            w = windingNumber(c, origB);
        }
        bool flip;
        if (shouldInclude(w, op, true, flip)) {
            DTri out = frag; if (flip) std::swap(out.v[1], out.v[2]);
            output.push_back(fromDTri(out, frag.mesh_id, false));
        }
    }
    for (const auto& frag : sr.splitB) {
        D3 c = d3_centroid(frag);
        int w = windingNumber(c, origA);
        bool flip;
        if (shouldInclude(w, op, false, flip)) {
            DTri out = frag; if (flip) std::swap(out.v[1], out.v[2]);
            output.push_back(fromDTri(out, frag.mesh_id, false));
        }
    }
    return output;
}

// ════════════════════════════════════════════════════════════════════════════
// §12  Pipeline stubs
// ════════════════════════════════════════════════════════════════════════════

void CherchiBackend::detectIntersections(const PolygonSoup&,std::vector<TriangleIntersection>&){}
void CherchiBackend::buildLocalArrangements(const PolygonSoup&,const std::vector<TriangleIntersection>&,std::vector<LocalArrangement>&){}
void CherchiBackend::classifyCells(const std::vector<LocalArrangement>&,BooleanOp,std::vector<bool>&){}
void CherchiBackend::extractResult(const PolygonSoup&,const std::vector<LocalArrangement>&,const std::vector<bool>&,BackendResult&){}

// ════════════════════════════════════════════════════════════════════════════
// §13  Main execution
// ════════════════════════════════════════════════════════════════════════════

BackendResult CherchiBackend::execute(const PolygonSoup& soup, BooleanOp op)
{
    BackendResult result;
    auto t0 = std::chrono::high_resolution_clock::now();

    if (soup.triangles.empty() || soup.mesh_tri_count.size() < 2) {
        result.success = false;
        result.error_message = "Cherchi: need two non-empty meshes";
        return result;
    }

    const uint32_t nA = soup.mesh_tri_count[0];
    const uint32_t nB = soup.mesh_tri_count[1];

    std::vector<DTri> origA(nA), origB(nB);
    for (uint32_t i = 0; i < nA; ++i) origA[i] = toDTri(soup.triangles[i]);
    for (uint32_t j = 0; j < nB; ++j) origB[j] = toDTri(soup.triangles[nA + j]);

    // Detect if B is open (plane, grid, etc.)
    bool bIsOpen = isMeshOpen(origB);
    FitPlane bPlane = {};
    if (bIsOpen) {
        bPlane = fitPlaneToMesh(origB);
    }

    // Split both meshes
    SplitResult sr = splitBothMeshes(soup, nA, nB, origA, origB);

    // Classify
    std::vector<Triangle> output;
    if (op == BooleanOp::Shatter) {
        output = classifyShatter(sr, origA, origB, bIsOpen, bPlane);
    } else {
        output = classifyStandard(sr, origA, origB, op, bIsOpen, bPlane);
    }

    // Write result
    result.output_soup.triangles = std::move(output);
    result.output_soup.quantization_scale[0] = soup.quantization_scale[0];
    result.output_soup.quantization_scale[1] = soup.quantization_scale[1];
    result.output_soup.quantization_scale[2] = soup.quantization_scale[2];
    result.output_soup.quantization_inv_scale[0] = soup.quantization_inv_scale[0];
    result.output_soup.quantization_inv_scale[1] = soup.quantization_inv_scale[1];
    result.output_soup.quantization_inv_scale[2] = soup.quantization_inv_scale[2];
    result.output_soup.per_axis_quantization = soup.per_axis_quantization;

    result.success = true;
    result.num_triangles = static_cast<uint32_t>(result.output_soup.triangles.size());

    auto t1 = std::chrono::high_resolution_clock::now();
    result.execution_time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return result;
}

} // namespace ember
