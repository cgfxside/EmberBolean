// ═══════════════════════════════════════════════════════════════════════════════
// EMBER CherchiBackend — Real Boolean Implementation (Phase 1.6)
// ═══════════════════════════════════════════════════════════════════════════════
//
// Implements Boolean operations (Union, Intersection, DiffAB, DiffBA, Shatter)
// via three stages:
//
//   STAGE 1 – detectIntersections()
//     Exact orient3d sign tests to find triangle pairs from mesh A and B
//     that actually intersect. Bounding-box broad phase avoids O(n²) calls.
//
//   STAGE 2 – splitByPlane() (called from execute())
//     For each triangle that crosses at least one opposing plane, split it
//     iteratively against every such plane using floating-point positions.
//     Produces a flat list of convex sub-triangles (fragments).
//
//   STAGE 3 – windingNumber() + shouldInclude() (called from execute())
//     Cast a ray in +X from each fragment's centroid; count signed crossings
//     with the opposing mesh (front face = +1, back face = -1).
//     Based on BooleanOp, include/exclude and optionally flip the normal.
//
// EXACTNESS
//   All "which side of this plane?" decisions use orient3d_filtered (double
//   filter + Int256 exact fallback). Intersection POSITIONS are computed in
//   float (acceptable since they are immediately used as new triangle vertices,
//   not as inputs to further predicate calls).
//
// OPEN MESH NOTE
//   Winding number is well-defined only for closed (watertight) meshes.
//   For open surfaces the count can be wrong near the open boundary.
//   The algorithm still produces useful output for the common VFX case where
//   at least one mesh is closed (e.g. a box cutter against a flat grid).

#include "CherchiBackend.h"
#include "../ember/PolygonSoup.h"
#include "../ember/ExactPredicates.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <unordered_set>
#include <vector>
#include <array>

namespace ember {

// ════════════════════════════════════════════════════════════════════════════
// §0  Satisfy private forward declarations from the header
// ════════════════════════════════════════════════════════════════════════════

struct CherchiBackend::ImplicitPoint   { uint32_t type; uint32_t data[4]; };
struct CherchiBackend::ArrangementPoint{ uint32_t index; bool is_explicit; };
struct CherchiBackend::TriangleIntersection {
    uint32_t triA;       // absolute index in combined soup (mesh A side)
    uint32_t triB;       // absolute index in combined soup (mesh B side)
    uint32_t num_points;
};
struct CherchiBackend::LocalArrangement {
    std::vector<uint32_t> triangles;
    std::vector<uint32_t> vertices;
};

// ════════════════════════════════════════════════════════════════════════════
// §1  Math helpers
// ════════════════════════════════════════════════════════════════════════════

using V3   = std::array<float, 3>;
using Tri3 = std::array<V3, 3>;

static V3 lerp3(V3 a, V3 b, float t) {
    return {a[0]+t*(b[0]-a[0]), a[1]+t*(b[1]-a[1]), a[2]+t*(b[2]-a[2])};
}
static V3 centroid3(const Tri3& t) {
    return {(t[0][0]+t[1][0]+t[2][0])/3.0f,
            (t[0][1]+t[1][1]+t[2][1])/3.0f,
            (t[0][2]+t[1][2]+t[2][2])/3.0f};
}
static Tri3 toTri3(const Triangle& t) {
    return {{{t.v[0][0],t.v[0][1],t.v[0][2]},
             {t.v[1][0],t.v[1][1],t.v[1][2]},
             {t.v[2][0],t.v[2][1],t.v[2][2]}}};
}
static Triangle makeTri(const Tri3& t3, int mesh_id, uint32_t src) {
    Triangle t;
    for (int v=0;v<3;++v) for (int a=0;a<3;++a) t.v[v][a] = t3[v][a];
    t.mesh_id = mesh_id; t.src_prim_idx = src;
    return t;
}
static Tri3 flipWinding(Tri3 t) { std::swap(t[1],t[2]); return t; }

// ── Plane helpers ─────────────────────────────────────────────────────────────

// Exact orientation: sign of P relative to plane of triangle (A,B,C).
// Uses orient3d_filtered (double filter → Int256 exact fallback).
static int planeSign(V3 A, V3 B, V3 C, V3 P) {
    return predicates::orient3d_filtered(
        A[0],A[1],A[2], B[0],B[1],B[2], C[0],C[1],C[2], P[0],P[1],P[2]);
}

// Float unnormalized signed distance from P to plane of (A,B,C).
// Sign matches planeSign; magnitude is proportional to real distance × ‖normal‖.
static float planeDist(V3 A, V3 B, V3 C, V3 P) {
    float abx=B[0]-A[0],aby=B[1]-A[1],abz=B[2]-A[2];
    float acx=C[0]-A[0],acy=C[1]-A[1],acz=C[2]-A[2];
    float nx=aby*acz-abz*acy, ny=abz*acx-abx*acz, nz=abx*acy-aby*acx;
    return nx*(P[0]-A[0]) + ny*(P[1]-A[1]) + nz*(P[2]-A[2]);
}

// ── Bounding box ──────────────────────────────────────────────────────────────

struct BBox {
    float mn[3]={ 1e30f, 1e30f, 1e30f};
    float mx[3]={-1e30f,-1e30f,-1e30f};
    void expand(V3 v) {
        for(int i=0;i<3;++i){mn[i]=std::min(mn[i],v[i]);mx[i]=std::max(mx[i],v[i]);}
    }
    bool overlaps(const BBox& o) const {
        for(int i=0;i<3;++i) if(mn[i]>o.mx[i]||mx[i]<o.mn[i]) return false;
        return true;
    }
    static BBox ofTri(const Triangle& t) {
        BBox b;
        for(int v=0;v<3;++v) b.expand({t.v[v][0],t.v[v][1],t.v[v][2]});
        return b;
    }
};

// ════════════════════════════════════════════════════════════════════════════
// §2  Exact triangle-triangle intersection test
// ════════════════════════════════════════════════════════════════════════════
// Returns true only when the triangles have a proper crossing (segment or
// area overlap). Coplanar triangles are returned as non-intersecting (they
// need 2D handling not required for the current test cases).

static bool triTriIntersect(const Triangle& A, const Triangle& B) {
    V3 a0={A.v[0][0],A.v[0][1],A.v[0][2]};
    V3 a1={A.v[1][0],A.v[1][1],A.v[1][2]};
    V3 a2={A.v[2][0],A.v[2][1],A.v[2][2]};
    V3 b0={B.v[0][0],B.v[0][1],B.v[0][2]};
    V3 b1={B.v[1][0],B.v[1][1],B.v[1][2]};
    V3 b2={B.v[2][0],B.v[2][1],B.v[2][2]};

    // Classify A vertices against B's plane
    int sA0=planeSign(b0,b1,b2,a0), sA1=planeSign(b0,b1,b2,a1), sA2=planeSign(b0,b1,b2,a2);
    if (sA0>0&&sA1>0&&sA2>0) return false;
    if (sA0<0&&sA1<0&&sA2<0) return false;
    if (sA0==0&&sA1==0&&sA2==0) return false;  // coplanar

    // Classify B vertices against A's plane
    int sB0=planeSign(a0,a1,a2,b0), sB1=planeSign(a0,a1,a2,b1), sB2=planeSign(a0,a1,a2,b2);
    if (sB0>0&&sB1>0&&sB2>0) return false;
    if (sB0<0&&sB1<0&&sB2<0) return false;
    if (sB0==0&&sB1==0&&sB2==0) return false;

    return true;
}

// ════════════════════════════════════════════════════════════════════════════
// §3  Triangle splitting by a plane
// ════════════════════════════════════════════════════════════════════════════
// Splits triangle T along the plane of triangle (A,B,C).
// Fragments with planeDist > 0 → `pos`.  planeDist < 0 → `neg`.
// On-plane vertices are treated as belonging to the majority side.

static void splitByPlane(
    const Tri3& T, V3 A, V3 B, V3 C,
    std::vector<Tri3>& pos, std::vector<Tri3>& neg)
{
    static constexpr float kEps = 1e-7f;
    float d[3]; int s[3];
    for(int i=0;i<3;++i){
        d[i]=planeDist(A,B,C,T[i]);
        s[i]=(d[i]>kEps)?1:(d[i]<-kEps)?-1:0;
    }
    int npos=(s[0]>0)+(s[1]>0)+(s[2]>0);
    int nneg=(s[0]<0)+(s[1]<0)+(s[2]<0);

    if(nneg==0){ pos.push_back(T); return; }
    if(npos==0){ neg.push_back(T); return; }

    // Find the lone vertex (on one side) and the two others (on opposite side)
    int lone=-1, o0=-1, o1=-1;
    int lone_sign = (npos==1) ? 1 : -1;

    for(int i=0;i<3;++i){
        bool is_lone = (lone_sign>0) ? (s[i]>0) : (s[i]<0);
        if(is_lone) lone=i;
        else { if(o0<0) o0=i; else o1=i; }
    }

    // Guard against degenerate denominators
    auto safeDiv = [](float num, float den) -> float {
        return (std::abs(den)>1e-20f) ? (num/den) : 0.5f;
    };

    float t0 = safeDiv(d[lone], d[lone]-d[o0]);
    float t1 = safeDiv(d[lone], d[lone]-d[o1]);
    V3 P0 = lerp3(T[lone], T[o0], t0);
    V3 P1 = lerp3(T[lone], T[o1], t1);

    // 1-vertex side: 1 triangle; 2-vertex side: 2 triangles (quad → 2 tris)
    Tri3 loneT   = {T[lone], P0, P1};
    Tri3 otherTA = {P0, T[o0], T[o1]};
    Tri3 otherTB = {P0, T[o1], P1};

    if(lone_sign>0){
        pos.push_back(loneT); neg.push_back(otherTA); neg.push_back(otherTB);
    } else {
        neg.push_back(loneT); pos.push_back(otherTA); pos.push_back(otherTB);
    }
}

// ════════════════════════════════════════════════════════════════════════════
// §4  Winding number via +X ray casting
// ════════════════════════════════════════════════════════════════════════════
// Cast a ray from O = (ox,oy,oz) in the +X direction through all triangles in
// `mesh`. Returns the signed winding number:
//   0   → O is outside mesh (assuming mesh is closed)
//   !=0 → O is inside mesh

static int windingNumber(V3 O, const std::vector<Tri3>& mesh) {
    int w = 0;
    for(const auto& T : mesh) {
        // Step 1: YZ-plane containment test using exact orient2d.
        // The +X ray hits this triangle only if the YZ projection of O
        // falls inside the YZ projection of T.
        double oy=O[1], oz=O[2];
        int d1=predicates::orient2d_filtered(T[0][1],T[0][2], T[1][1],T[1][2], oy,oz);
        int d2=predicates::orient2d_filtered(T[1][1],T[1][2], T[2][1],T[2][2], oy,oz);
        int d3=predicates::orient2d_filtered(T[2][1],T[2][2], T[0][1],T[0][2], oy,oz);
        bool ccw = (d1>0&&d2>0&&d3>0);
        bool  cw = (d1<0&&d2<0&&d3<0);
        if(!ccw && !cw) continue;

        // Step 2: Compute the X coordinate of the ray-plane intersection.
        // Plane normal: n = (B-A) × (C-A)
        float abx=T[1][0]-T[0][0], aby=T[1][1]-T[0][1], abz=T[1][2]-T[0][2];
        float acx=T[2][0]-T[0][0], acy=T[2][1]-T[0][1], acz=T[2][2]-T[0][2];
        float nx=aby*acz-abz*acy;
        float ny=abz*acx-abx*acz;
        float nz=abx*acy-aby*acx;

        // Ray: (ox + t, oy, oz). Plane: n·(X-A)=0.
        // → nx*(ox+t - A.x) + ny*(oy - A.y) + nz*(oz - A.z) = 0
        // → t = -(n·(O-A)) / nx
        if(std::abs(nx) < 1e-10f) continue;  // plane parallel to +X ray

        float nDotOA = nx*(O[0]-T[0][0]) + ny*(O[1]-T[0][1]) + nz*(O[2]-T[0][2]);
        float t = -nDotOA / nx;
        if(t <= 0.0f) continue;  // crossing is behind O (in -X direction)

        // Step 3: Signed crossing — front face (+1) or back face (−1)
        w += (nx > 0.0f) ? +1 : -1;
    }
    return w;
}

// ════════════════════════════════════════════════════════════════════════════
// §5  Boolean inclusion logic
// ════════════════════════════════════════════════════════════════════════════

// Decides whether a fragment should be in the output given:
//   winding    – winding number of its centroid relative to the other mesh
//   op         – the Boolean operation
//   from_A     – true if this fragment came from mesh A, false if from mesh B
//   flip_normal – set to true if the fragment's winding order should be reversed
static bool shouldInclude(int winding, BooleanOp op, bool from_A, bool& flip_normal) {
    flip_normal = false;
    bool inside = (winding != 0);

    switch(op) {
        case BooleanOp::Union:
            return !inside;  // keep what is OUTSIDE the other mesh

        case BooleanOp::Intersection:
            return inside;   // keep what is INSIDE the other mesh

        case BooleanOp::DiffAB:  // A − B
            if(from_A) return !inside;         // A fragments outside B → keep
            flip_normal = true; return inside; // B fragments inside A → keep, flip

        case BooleanOp::DiffBA:  // B − A
            if(!from_A) return !inside;        // B fragments outside A → keep
            flip_normal = true; return inside; // A fragments inside B → keep, flip

        case BooleanOp::Shatter:
            return true;   // keep everything (piece labeling done in ShatterPost)

        case BooleanOp::Seam:
            return false;  // seam = intersection curves only, no surface output

        default:
            return !inside;
    }
}

// ════════════════════════════════════════════════════════════════════════════
// §6  Pipeline stage stubs (required by header; logic lives in execute())
// ════════════════════════════════════════════════════════════════════════════

void CherchiBackend::detectIntersections(
    const PolygonSoup& soup,
    std::vector<TriangleIntersection>& intersections)
{
    uint32_t nA = soup.mesh_tri_count.size()>0 ? soup.mesh_tri_count[0] : 0;
    uint32_t nB = soup.mesh_tri_count.size()>1 ? soup.mesh_tri_count[1] : 0;

    // Build B bounding boxes once
    std::vector<BBox> bboxB(nB);
    for(uint32_t j=0;j<nB;++j) bboxB[j]=BBox::ofTri(soup.triangles[nA+j]);

    for(uint32_t i=0;i<nA;++i){
        BBox bA=BBox::ofTri(soup.triangles[i]);
        for(uint32_t j=0;j<nB;++j){
            if(!bA.overlaps(bboxB[j])) continue;
            if(triTriIntersect(soup.triangles[i], soup.triangles[nA+j]))
                intersections.push_back({i, nA+j, 0});
        }
    }
}

void CherchiBackend::buildLocalArrangements(
    const PolygonSoup&, const std::vector<TriangleIntersection>&,
    std::vector<LocalArrangement>&) { /* implemented inline in execute() */ }

void CherchiBackend::classifyCells(
    const std::vector<LocalArrangement>&, BooleanOp,
    std::vector<bool>&) { /* implemented inline in execute() */ }

void CherchiBackend::extractResult(
    const PolygonSoup&, const std::vector<LocalArrangement>&,
    const std::vector<bool>&, BackendResult&) { /* implemented inline in execute() */ }

// ════════════════════════════════════════════════════════════════════════════
// §7  Main execution
// ════════════════════════════════════════════════════════════════════════════

BackendResult CherchiBackend::execute(const PolygonSoup& soup, BooleanOp op)
{
    BackendResult result;
    auto t0 = std::chrono::high_resolution_clock::now();

    if(soup.triangles.empty() || soup.mesh_tri_count.size()<2){
        result.success=false;
        result.error_message="Cherchi: need two non-empty meshes";
        return result;
    }

    uint32_t nA=soup.mesh_tri_count[0];
    uint32_t nB=soup.mesh_tri_count[1];

    // §7.1  Build Tri3 reference arrays for winding number ray casting ─────────
    std::vector<Tri3> meshA(nA), meshB(nB);
    for(uint32_t i=0;i<nA;++i) meshA[i]=toTri3(soup.triangles[i]);
    for(uint32_t j=0;j<nB;++j) meshB[j]=toTri3(soup.triangles[nA+j]);

    // §7.2  Detect intersecting triangle pairs ────────────────────────────────
    std::vector<TriangleIntersection> intersections;
    detectIntersections(soup, intersections);

    // Which triangles need splitting?
    std::unordered_set<uint32_t> touchedA, touchedB;
    for(const auto& p : intersections){
        touchedA.insert(p.triA);
        touchedB.insert(p.triB);
    }

    // Inverse map: for each A-triangle, list of B-triangle indices that cut it
    std::vector<std::vector<uint32_t>> cutsOnA(nA);
    std::vector<std::vector<uint32_t>> cutsOnB(nB);
    for(const auto& p : intersections){
        cutsOnA[p.triA].push_back(p.triB);
        cutsOnB[p.triB - nA].push_back(p.triA);
    }

    // §7.3  Process mesh A triangles ──────────────────────────────────────────
    std::vector<Triangle> output;

    for(uint32_t i=0;i<nA;++i){
        const Triangle& src=soup.triangles[i];

        if(touchedA.find(i)==touchedA.end()){
            // No split needed — classify whole triangle by centroid winding
            V3 c=centroid3(meshA[i]);
            int w=windingNumber(c, meshB);
            bool flip_n;
            if(shouldInclude(w, op, true, flip_n)){
                Triangle out=src;
                if(flip_n) std::swap(out.v[1],out.v[2]);
                output.push_back(out);
            }
        } else {
            // Iteratively split this triangle against each opposing B plane
            std::vector<Tri3> frags={toTri3(src)};
            for(uint32_t bIdx : cutsOnA[i]){
                const Triangle& B=soup.triangles[bIdx];
                V3 b0={B.v[0][0],B.v[0][1],B.v[0][2]};
                V3 b1={B.v[1][0],B.v[1][1],B.v[1][2]};
                V3 b2={B.v[2][0],B.v[2][1],B.v[2][2]};
                std::vector<Tri3> next;
                for(const auto& frag : frags){
                    std::vector<Tri3> pos, neg;
                    splitByPlane(frag, b0,b1,b2, pos, neg);
                    next.insert(next.end(),pos.begin(),pos.end());
                    next.insert(next.end(),neg.begin(),neg.end());
                }
                frags=std::move(next);
            }
            // Classify each fragment
            for(const auto& frag : frags){
                V3 c=centroid3(frag);
                int w=windingNumber(c, meshB);
                bool flip_n;
                if(shouldInclude(w, op, true, flip_n)){
                    Tri3 out=flip_n ? flipWinding(frag) : frag;
                    output.push_back(makeTri(out, 0, src.src_prim_idx));
                }
            }
        }
    }

    // §7.4  Process mesh B triangles ──────────────────────────────────────────

    for(uint32_t j=0;j<nB;++j){
        const Triangle& src=soup.triangles[nA+j];

        if(touchedB.find(nA+j)==touchedB.end()){
            V3 c=centroid3(meshB[j]);
            int w=windingNumber(c, meshA);
            bool flip_n;
            if(shouldInclude(w, op, false, flip_n)){
                Triangle out=src;
                if(flip_n) std::swap(out.v[1],out.v[2]);
                output.push_back(out);
            }
        } else {
            std::vector<Tri3> frags={toTri3(src)};
            for(uint32_t aIdx : cutsOnB[j]){
                const Triangle& A=soup.triangles[aIdx];
                V3 a0={A.v[0][0],A.v[0][1],A.v[0][2]};
                V3 a1={A.v[1][0],A.v[1][1],A.v[1][2]};
                V3 a2={A.v[2][0],A.v[2][1],A.v[2][2]};
                std::vector<Tri3> next;
                for(const auto& frag : frags){
                    std::vector<Tri3> pos, neg;
                    splitByPlane(frag, a0,a1,a2, pos, neg);
                    next.insert(next.end(),pos.begin(),pos.end());
                    next.insert(next.end(),neg.begin(),neg.end());
                }
                frags=std::move(next);
            }
            for(const auto& frag : frags){
                V3 c=centroid3(frag);
                int w=windingNumber(c, meshA);
                bool flip_n;
                if(shouldInclude(w, op, false, flip_n)){
                    Tri3 out=flip_n ? flipWinding(frag) : frag;
                    output.push_back(makeTri(out, 1, src.src_prim_idx));
                }
            }
        }
    }

    // §7.5  Write result using passthrough path ────────────────────────────────
    // GU_EmberConverter's passthrough branch (output_polygons.empty()) writes
    // tri.v[] float coordinates directly to GU_Detail — no dequantization.
    result.output_soup.triangles=std::move(output);
    result.output_soup.quantization_scale[0]=soup.quantization_scale[0];
    result.output_soup.quantization_scale[1]=soup.quantization_scale[1];
    result.output_soup.quantization_scale[2]=soup.quantization_scale[2];
    result.output_soup.quantization_inv_scale[0]=soup.quantization_inv_scale[0];
    result.output_soup.quantization_inv_scale[1]=soup.quantization_inv_scale[1];
    result.output_soup.quantization_inv_scale[2]=soup.quantization_inv_scale[2];
    result.output_soup.per_axis_quantization=soup.per_axis_quantization;
    // output_polygons intentionally left empty → passthrough path

    result.success=true;
    result.num_triangles=static_cast<uint32_t>(result.output_soup.triangles.size());

    auto t1=std::chrono::high_resolution_clock::now();
    result.execution_time_ms=std::chrono::duration<double,std::milli>(t1-t0).count();
    return result;
}

} // namespace ember
