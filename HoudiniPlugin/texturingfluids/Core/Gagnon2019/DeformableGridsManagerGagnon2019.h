#ifndef __DeformableGridGagnon2019_h__
#define __DeformableGridGagnon2019_h__

#include <Math/Vec3.h>
#include <Core/Gagnon2019/ParticleTrackerManagerGagnon2019.h>
#include "Core/Deformations/DistortionMetricSorkine2002.h"
#include <GEO/GEO_PointTree.h>
#include <GU/GU_RayIntersect.h>

namespace TexturingFluids {

#define VERBOSE 0


class DeformableGridsManagerGagnon2019 : public ParticleTrackerManagerGagnon2019
{
public:

    DeformableGridsManagerGagnon2019(GU_Detail *surfaceGdp, GU_Detail *trackersGdp);


    void CreateGridsBasedOnMesh(GU_Detail *gdp,GU_Detail *surfaceGdp, GU_Detail *trackersGdp, ParametersDeformablePatches params,vector<GA_Offset> trackers,  GEO_PointTreeGAOffset &tree);
    void CreateGridBasedOnMesh(GU_Detail *deformableGridsGdp,GU_Detail *surfaceGdp, GU_Detail *trackersGdp, ParametersDeformablePatches params, GA_Offset ppt,  GEO_PointTreeGAOffset &tree);

    void AdvectGrids(GU_Detail *gdp, GU_Detail *trackersGdp, ParametersDeformablePatches params, GEO_PointTreeGAOffset &tree, GU_Detail *surfaceGdp);

    void ConnectivityTest(const GU_Detail *gdp,GA_Offset point, GA_PointGroup *grp, set<GA_Offset> &pointsAround,set<GA_Offset> &group);

    bool UVFlattening(GU_Detail &tempGdp, GU_Detail *trackersGdp, GU_Detail *deformableGridsGdp,
                      GA_Offset tracker, GA_Offset closestPoint,
                      GA_PointGroup *pointGroup, GA_PointGroup *tempPointGroup,
                      set<GA_Offset> &pointsAround,
                      float scaling,
                      ParametersDeformablePatches params);
    void FlagBoundaries(GU_Detail *deformableGridsGdp);

    const string gridGroupName = "grids";

    double  uvFlatteningTime;
    int     nbOfFlattenedPatch;
    double  gridMeshCreation;
    double  gridAdvectionTime;
    const string approachName   = "[Deformable Grids Gagno 2019]";

    int numberOfDegeneratedGrid;


protected :

    //void ConnectivityTest(const GU_Detail *gdp,GA_Offset point, GA_PointGroup *grp, GEO_PointTreeGAOffset::IdxArrayType pointsAround,set<GA_Offset> &group);


    const string uvArrayName = "uvs";
    const string alphaArrayName = "alphas";
    const string patchIdsName = "patchIds";

    const char*    alpha0Name      = "alpha0";
    const char*    distortionWeightName  = "w";
    const char*    vertexWeightName  = "vw";
    const char*    distortionWeight0Name  = "w0";
    const char*    disrortionMinThreshold = "t_min";
    const char*    distortionMaxThreshold = "t_max";
    const char*    initialVertexAngle = "v_a0";
    const char*    uvName = "uvw";

    const string triangleArrayName = "triangleRef";
    GA_Attribute        *triangleRefAtt;
    const GA_AIFNumericArray *triangleArray;

    vector<set<map<GA_Offset, UT_Vector3> > > uvs;
    map<int,UT_Vector3> gridCenterPosition;
    bool useUvFlattening = true;



};

}

#endif
