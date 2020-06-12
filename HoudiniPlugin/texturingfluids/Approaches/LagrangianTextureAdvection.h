#ifndef __Yu2011ExtentionInterface_h__
#define __Yu2011ExtentionInterface_h__
#include <cassert>
#include <cmath>
#include <vector>
#include <SOP/SOP_Node.h>
#include <GEO/GEO_PrimPart.h>
#include <Math/Vec3.h>
#include "Images/Image.h"
#include "Set/SpatialGrid.h"
#include <GU/GU_Flatten.h>
#include <Strategies/StrategyPatchSurfaceSynthesis.h>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/tbb.h>

namespace TexturingFluids {

class LagrangianTextureAdvection
{

public:

    //=========================== BUILD ======================================
    LagrangianTextureAdvection();
    ~LagrangianTextureAdvection();

    //==========================================================================
    void Synthesis(GU_Detail* gdp, GU_Detail* surface, GU_Detail* trackersGdp, GU_Detail* levelSet, GU_Detail *surfaceLowRes, ParametersDeformablePatches params);

private :

	
};
}

#endif
