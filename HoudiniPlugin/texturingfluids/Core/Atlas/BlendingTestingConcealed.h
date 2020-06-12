#ifndef __BlendingTestingConcealed_h__
#define __BlendingTestingConcealed_h__

#include <string>
#include <vector>
#include "Blending.h"


namespace TexturingFluids {

class BlendingTestingConcealed : public Blending
{

public:

    static Pixel Blend(GU_Detail* deformableGrids, int i, int j, float w, float h,
                       int pixelPositionX, int pixelPositionY,
                       vector<int> &sortedPatches,
                       vector<UT_Vector3> &surfaceUv,
                       vector<UT_Vector3> &surfacePosition,
                       map<int,UT_Vector3> &trackersNormal,
                       map<int,UT_Vector3> &trackersPosition,
                       map<int,UT_Vector3> &trackersUVPosition,
                       map<int, bool> &usePatches,
                       map<string,GU_RayIntersect*> &rays,
                       GA_RWHandleV3 &attPointUV,
                       map<int,float> &patchBlend,
                       vector<ImageCV*> textureExemplars,

                       ParametersDeformablePatches params);


};
}

#endif
