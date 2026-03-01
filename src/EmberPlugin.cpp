/**
 * @file EmberPlugin.cpp
 * @brief DSO entry point for EMBER Boolean Houdini plugin
 */

#include "SOP_EmberBoolean.h"

#include <UT/UT_DSOVersion.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>

void newSopOperator(OP_OperatorTable* table)
{
    table->addOperator(new OP_Operator(
        "ember_boolean",                        // internal name
        "EMBER Boolean",                        // UI label
        SOP_EmberBoolean::myConstructor,        // constructor
        SOP_EmberBoolean::buildTemplates(),     // parameters (public accessor)
        2,                                      // min inputs
        2,                                      // max inputs
        nullptr,                                // variables
        OP_FLAG_GENERATOR                       // flags
    ));
}
