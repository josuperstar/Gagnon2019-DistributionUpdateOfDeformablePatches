#include "LagrangianTextureAdvection.h"

#include <fstream>
#include <vector>
#include <algorithm>
#include <SYS/SYS_Math.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Matrix4.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <PRM/PRM_Include.h>
#include <PRM/PRM_SpareData.h>
#include <SOP/SOP_Guide.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <GA/GA_ElementWrangler.h>
#include <algorithm>
#include <ctime>
#include <GEO/GEO_PointTree.h>
#include <GU/GU_NeighbourList.h>

#include <GU/GU_Flatten.h>
#include <GU/GU_RayIntersect.h>

#include <Core/Gagnon2019/PatchedSurfaceGagnon2019.h>
#include <Core/Gagnon2019/Bridson2012PoissonDiskDistribution.h>

LagrangianTextureAdvection::LagrangianTextureAdvection()
{
}

LagrangianTextureAdvection::~LagrangianTextureAdvection()
{
}

void LagrangianTextureAdvection::Synthesis(GU_Detail *gdp, GU_Detail *surfaceGdp, GU_Detail *trackersGdp, GU_Detail *levelSet, GU_Detail *surfaceLowResGdp,  ParametersDeformablePatches params)
{
    PatchedSurfaceGagnon2019 surface(surfaceGdp, trackersGdp);
    cout << "[LagrangianTextureAdvection::Synthesis] "<<params.frame<<endl;
    //params.useDynamicTau = false;

    std::clock_t start;
    start = std::clock();
    vector<GA_Offset> newPatchesPoints;
    vector<GA_Offset> trackers;
    cout << "reference gdp created"<<endl;
    const GA_SaveOptions *options;
    UT_StringArray *errors;

    GA_PointGroup *surfaceGroup = (GA_PointGroup *)surfaceGdp->pointGroups().find(surface.surfaceGroupName.c_str());
    if (surfaceGroup == 0x0)
    {
        cout << "There is no surface group to synthesis"<<endl;
        return;
    }
    //=======================================================
    GA_PointGroup *grp = (GA_PointGroup *)gdp->pointGroups().find(surface.markerGroupName.c_str());

    GU_RayIntersect ray(gdp);
    ray.init();
    GEO_PointTreeGAOffset surfaceTree;
    surfaceTree.build(surfaceGdp, NULL);
    GEO_PointTreeGAOffset surfaceLowResTree;
    surfaceLowResTree.build(surfaceLowResGdp, NULL);

    //=========================== CORE ALGORITHM ============================


    //---- for visualisation purpose

    //string beforeUpdateString = params.trackersFilename + "beforeAdvection.bgeo";
    //const char* filename = beforeUpdateString.c_str();//"dlttest.bgeo";
    //trackersGdp->save(filename,options,errors);
    //----------------------------------


    bool usingOnlyPoissonDisk = false;


    if(params.startFrame == params.frame)
    {
        surface.PoissonDiskSampling(levelSet,trackersGdp,params);
        surface.CreateAndUpdateTrackersBasedOnPoissonDisk(surfaceGdp,trackersGdp, surfaceGroup,params);
        if (!usingOnlyPoissonDisk)
            surface.CreateGridsBasedOnMesh(gdp,surfaceLowResGdp,trackersGdp, params,newPatchesPoints,surfaceLowResTree);
    }
    else
    {
        surface.AdvectSingleTrackers(surfaceLowResGdp,trackersGdp, params);
        if (!usingOnlyPoissonDisk)
            surface.AdvectGrids(gdp,trackersGdp,params,surfaceLowResTree,surfaceLowResGdp);
        if (params.updateDistribution)
        {
            surface.PoissonDiskSampling(levelSet,trackersGdp,params); //Poisson disk on the level set
            //surface.CreateAndUpdateTrackersBasedOnPoissonDisk(surfaceGdp,trackersGdp, surfaceGroup,params);
        }
        surface.CreateAndUpdateTrackersBasedOnPoissonDisk(surfaceGdp,trackersGdp, surfaceGroup,params);
        if (!usingOnlyPoissonDisk)
            surface.CreateGridsBasedOnMesh(gdp,surfaceLowResGdp,trackersGdp, params,newPatchesPoints,surfaceLowResTree);
        surface.DeleteUnusedPatches(gdp, trackersGdp,params);

    }
    if (!usingOnlyPoissonDisk)
    {
        //For the blending computation, we create uv array per vertex that we called patch
        surface.AddDeformablePatchesUsingBarycentricCoordinates(gdp, surfaceGdp,trackersGdp, params,surfaceTree,ray);
    }

    //=======================================================================

    cout << surface.approachName<<" Done"<<endl;
    cout << "Clear surface tree"<<endl;
    surfaceTree.clear();
    ray.clear();

    cout << surface.approachName<< " saving grids data"<<endl;
    const char* filenameGrids = params.deformableGridsFilename.c_str();//"dlttest.bgeo";
    gdp->save(filenameGrids,options,errors);

    cout << surface.approachName<< " saving trackers data"<<endl;
    const char* filenameTrackers = params.trackersFilename.c_str();//"dlttest.bgeo";
    trackersGdp->save(filenameTrackers,options,errors);

    //================================================================
    std::clock_t cleaningStart;
    cleaningStart = std::clock();
    cout<< "Clear, Destroy and merge"<<endl;
    gdp->clearAndDestroy();
    gdp->copy(*surfaceGdp);

    int nbPatches = surface.GetNumberOfPatches();

    float cleaningSurface = (std::clock() - cleaningStart) / (double) CLOCKS_PER_SEC;
    cout << "--------------------------------------------------------------------------------"<<endl;
    cout << surface.approachName<<" Poisson Disk Sampling "<<surface.poissondisk<<endl;
    cout << surface.approachName<<" Grid mesh on time "<<surface.gridMeshCreation<<endl;
    cout << surface.approachName<<" Uv flattening time "<<surface.uvFlatteningTime<<" for "<<surface.nbOfFlattenedPatch<<" patches"<<endl;
    cout << surface.approachName<<" Tracker advection time "<<surface.markerAdvectionTime<<endl;
    cout << surface.approachName<<" Grid advection time "<<surface.gridAdvectionTime<<endl;
    cout << surface.approachName<<" Patch creation time "<<surface.patchCreationTime<<endl;
    cout << surface.approachName<<" Clear and Destroy "<<cleaningSurface<<endl;
    cout << surface.approachName<<" Update distribution "<<surface.updatePatchesTime<<endl;

    float total = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    cout << surface.approachName<< " TOTAL: "<<total<<endl;

    std::ofstream outfile;
    outfile.open("core.csv", std::ios_base::app);
    outfile <<surface.poissondisk<<","<< surface.gridMeshCreation << ","<<surface.uvFlatteningTime << ","<<surface.markerAdvectionTime
            <<","<<surface.gridAdvectionTime<<","<<surface.patchCreationTime << ","<<surface.updatePatchesTime<<","<<nbPatches<<endl;

    cout << "--------------------------------------------------------------------------------"<<endl;
}

