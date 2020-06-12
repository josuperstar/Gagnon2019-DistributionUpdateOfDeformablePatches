#ifndef __AtlasInterface_h__
#define __AtlasInterface_h__

#include <cassert>
#include <cmath>
#include <vector>
#include <SOP/SOP_Node.h>
#include <GEO/GEO_PrimPart.h>
//#include <GA_ElementGroup.h>
#include <Math/Vec3.h>
#include "Images/Image.h"
#include "Set/SpatialGrid.h"


#include <GU/GU_Flatten.h>
#include <Core/Deformations/ParametersDeformablePatches.h>
#include <Strategies/StrategyPatchSurfaceSynthesis.h>

namespace TexturingFluids {

class AtlasInterface
{

public:

    //=========================== BUILD ======================================
    AtlasInterface();
    ~AtlasInterface();

    //==========================================================================

    bool Synthesis(GU_Detail* gdp,  GU_Detail *surfaceGdp, GU_Detail *trackersGdp, ParametersDeformablePatches params);


private :



};
}

#endif
