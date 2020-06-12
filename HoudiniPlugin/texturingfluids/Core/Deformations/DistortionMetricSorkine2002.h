#ifndef __YU2011DISTORTION_h__
#define __YU2011DISTORTION_h__
#include <cassert>
#include <cmath>
#include <vector>
#include <SOP/SOP_Node.h>
#include <GEO/GEO_PrimPart.h>
#include <Math/Vec3.h>
#include "Images/Image.h"
#include "Set/SpatialGrid.h"

#include <GU/GU_Flatten.h>
#include "ParametersDistortion.h"

namespace TexturingFluids {


class DistortionMetricSorkine2002
{

public:

   bool ComputeDistortion(GU_Detail *trackers, GU_Detail *gridGdp, GA_Offset trackerPpt, GA_PointGroup* pointGrp, GA_PrimitiveGroup *primGroup, ParametersDistortion params);

};
}

#endif
