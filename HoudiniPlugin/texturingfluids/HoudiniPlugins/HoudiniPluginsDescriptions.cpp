/**
 * Code de François Dagenais
 *
 */

#include <OP/OP_Director.h>


#include "Yu2011Plugin.h"
#include "AtlasYu2011Plugin.h"

// -----------------------------------------------------------------------------
// Add our plugins to Houdini's plugins list
// -----------------------------------------------------------------------------
void newSopOperator(OP_OperatorTable *table)
{
     table->addOperator(new OP_Operator("hdk_Yu2011",
                                        "LagrangianTextureAdvection",
                                        LagrangianTextureAdvectionPlugin::myConstructor,
                                        LagrangianTextureAdvectionPlugin::myTemplateList,
                                        5,
                                        5,
                                        0));


     table->addOperator(new OP_Operator("hdk_AtlasYu2011",
                                       "AtlasYu2011",
                                       AtlasYu2011Plugin::myConstructor,
                                       AtlasYu2011Plugin::myTemplateList,
                                       3,
                                       3,
                                       0));


}
