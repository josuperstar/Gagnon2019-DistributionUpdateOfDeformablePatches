#ifndef __PatchedSurface__
#define __PatchedSurface__

#include <Math/Vec3.h>
#include <Core/Gagnon2019/DeformableGridsManagerGagnon2019.h>
#include <GEO/GEO_PointTree.h>
#include <GU/GU_RayIntersect.h>

namespace TexturingFluids {

#define VERBOSE 0


class PatchedSurfaceGagnon2019 : public DeformableGridsManagerGagnon2019
{
public:

    PatchedSurfaceGagnon2019(GU_Detail *surface, GU_Detail *trackersGdp);
    ~PatchedSurfaceGagnon2019();

    void PoissonDiskSampling(GU_Detail *surfaceGdp, GU_Detail *trackers, ParametersDeformablePatches params);
    void AddDeformablePatcheUsingBarycentricCoordinates(GU_Detail *deformableGridsGdp,GU_Detail *surfaceGdp, GU_Detail *trackersGdp, GA_Offset ppt, ParametersDeformablePatches params, GEO_PointTreeGAOffset &surfaceTree,  GU_RayIntersect &ray);
    void AddDeformablePatchesUsingBarycentricCoordinates(GU_Detail *gdp, GU_Detail* surface, GU_Detail *trackersGdp, ParametersDeformablePatches params,  GEO_PointTreeGAOffset &surfaceTree, GU_RayIntersect &ray);
    void DeleteUnusedPatches(GU_Detail *gdp, GU_Detail *trackersGdp, ParametersDeformablePatches params);

    //for test purpose
    GA_Offset CreateAPatch(GU_Detail *trackers, UT_Vector3 position, UT_Vector3 normal, ParametersDeformablePatches params);

    double poissondisk;
    double  patchCreationTime;
    double  updatePatchesTime;
    const string approachName   = "[Patched Surface]";

private :

    //-------------- VERTEX ARRAY -----------------
    GA_Attribute        *uvsAtt;
    const GA_AIFNumericArray *uvsArray;

    GA_Attribute        *patchIdsArrayAttrib;
    const GA_AIFNumericArray *patchIdsAtt;

    GA_Attribute        *alphaArrayAtt;
    const GA_AIFNumericArray *alphaAtt;
    //---------------------------------------------

    map<string,GU_RayIntersect*> rays;


};

}

#endif
