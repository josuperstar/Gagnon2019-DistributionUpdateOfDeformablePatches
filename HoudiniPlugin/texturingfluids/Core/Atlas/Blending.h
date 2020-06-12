#ifndef __Blending_h__
#define __Blending_h__

#include <string>
#include <vector>
#include "Images/Color.h"
#include "Images/ImageCV.h"

#include <SOP/SOP_Node.h>
#include <GEO/GEO_PrimPart.h>
#include <GEO/GEO_PointTree.h>
#include <GU/GU_RayIntersect.h>
#include "Core/Deformations/ParametersDeformablePatches.h"


namespace TexturingFluids {

class Blending
{

public:

    static Pixel Blend(GU_Detail* trackersGdp, GU_Detail* deformableGrids, int i, int j, float w, float h,
                        int pixelPositionX, int pixelPositionY,
                        vector<int> &sortedPatches,
                        vector<UT_Vector3> &surfaceUv,
                        vector<UT_Vector3> &surfacePosition,
                        map<int,UT_Vector3> &trackersPosition,
                        bool useDeformableGrids,
                        map<string,GU_RayIntersect*> &rays,
                        map<int,Pixel> &patchColors,
                        Pixel alphaColor,
                        Pixel RM,           //Mean Value
                        GA_ROHandleF &attAlpha,
                        GA_RWHandleV3 &attPointUV,
                        int life,
                        map<int,float> &patchBlend,
                        map<int, vector<UT_Vector3> > &patchUvs,
                        map<int, vector<float> > &alphasMap,
                        vector<ImageCV*> textureExemplars,
                        ImageCV *textureExemplar1ImageMask,
                        ImageCV *displacementMapImage,
                        bool computeDisplacement,
                        bool renderColoredPatches,
                        Pixel &R1,
                        Pixel &displacementSum,
                        ParametersDeformablePatches params);


    static void Clamp(Pixel &pixel)
    {

        /*
        pixel.A = 1;
        pixel.R /= 2;
        pixel.G /= 2;
        pixel.B /= 2;


        pixel.R += 0.25;
        pixel.G += 0.25;
        pixel.B += 0.25;

        pixel.R -= 0.1;
        pixel.G -= 0.1;
        pixel.B -= 0.1;
        */

        //cout << pixel.A << " "<<pixel.R<<" "<<pixel.G << " "<<pixel.B<<endl;

        if (pixel.A > 1)
            pixel.A = 1;

        if (pixel.R > 1)
            pixel.R = 1;
        if (pixel.G > 1)
            pixel.G = 1;
        if (pixel.B > 1)
            pixel.B = 1;

        if (pixel.A < 0)
            pixel.A = 0;

        if (pixel.R < 0)
            pixel.R = 0;
        if (pixel.G < 0)
            pixel.G = 0;
        if (pixel.B < 0)
            pixel.B = 0;

    }

};
}

#endif
