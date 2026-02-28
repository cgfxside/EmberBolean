// ============================================================================
// EmberPlugin.cpp
// DSO Entry Point for EMBER Boolean Engine HDK Integration
// ============================================================================

#include <UT/UT_DSOVersion.h>
#include <OP/OP_OperatorTable.h>
#include <OP/OP_Operator.h>
#include <SOP/SOP_Node.h>

#include "SOP_EmberBoolean.h"

#define EMBER_PLUGIN_NAME    "ember_boolean"
#define EMBER_PLUGIN_LABEL   "Ember Boolean"

/// @brief Register the EMBER Boolean SOP with Houdini
void newSopOperator(OP_OperatorTable* table)
{
    if (!table) return;

    OP_Operator* op = new OP_Operator(
        EMBER_PLUGIN_NAME,
        EMBER_PLUGIN_LABEL,
        SOP_EmberBoolean::myConstructor,
        SOP_EmberBoolean::buildTemplates(),
        2,          // min inputs
        2,          // max inputs
        nullptr,    // local variables
        0           // flags (not a generator)
    );

    if (op) {
        // FIX: setOpVersion() and setOpDescription() do not exist in H20/H21.
        // Operator metadata is expressed only through the constructor arguments
        // and the PRM_Template list.
        table->addOperator(op);
    }
}
