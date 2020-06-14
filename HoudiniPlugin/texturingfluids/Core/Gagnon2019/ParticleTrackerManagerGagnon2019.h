#ifndef __ParticleTrackerGagnon2019_h__
#define __ParticleTrackerGagnon2019_h__

#include <Math/Vec3.h>
#include <Strategies/StrategyPatchSurfaceSynthesis.h>
#include <GEO/GEO_PointTree.h>
#include <GU/GU_RayIntersect.h>

namespace TexturingFluids {

#define VERBOSE 0


class ParticleTrackerManagerGagnon2019
{
public:

    ParticleTrackerManagerGagnon2019(GU_Detail *surfaceGdp, GU_Detail *trackersGdp);
    void CreateAndUpdateTrackersBasedOnPoissonDisk(GU_Detail* surface,GU_Detail *trackers, GA_PointGroup *surfaceGroup, ParametersDeformablePatches params);
    void CreateAndUpdateTrackerBasedOnPoissonDisk(GU_Detail *surface, GU_Detail *trackersGdp, GA_Offset ppt, GA_PointGroup *surfaceGroup,  ParametersDeformablePatches params);
    void UpdateTrackersAndTangeant(GU_Detail* surface,GU_Detail *trackers, GA_PointGroup *surfaceGroup, ParametersDeformablePatches params);
    void AdvectSingleTrackers(GU_Detail *surfaceGdp, GU_Detail *trackers, ParametersDeformablePatches params);
    void AdvectTrackersAndTangeants(GU_Detail *surfaceGdp, GU_Detail *trackers, ParametersDeformablePatches params);
    void ComputeDivergence(GU_Detail *surfaceGdp, GU_Detail *trackers, ParametersDeformablePatches params, GEO_PointTreeGAOffset &tree);
    void ComputeDensity(GU_Detail *surfaceGdp, GU_Detail *trackers, ParametersDeformablePatches params, GEO_PointTreeGAOffset &tree);
    void DeleteTracker(GU_Detail* trackers,int trackerId);

    int GetNumberOfPatches(){return numberOfPatches;}

    int NumberOfPatchesToDelete(GU_Detail *trackersGdp);

    const string markerGroupName = "markers";
    const string surfaceGroupName = "surface";
    const string approachName   = "[Particle Tracker Gagnon 2019]";

    double  markerAdvectionTime;

    int numberOfPatches;
    int numberOfInitialPatchFlagToDelete;
    int numberOfConcealedPatches;
    int numberOfNewPatches;
    int numberOfDetachedPatches;
    int numberOfLonelyTracker;
    int numberOfNewAndLonelyTracker;
    int numberOfInitialPatches;
    int numberOfDistortedPatches;


protected :

    bool tackerPolygon = false;

    int maxId = 0;
    const char*    randomThresholdDistortion = "t_rand";
    float epsilon = 0.0001;

    GA_RWHandleV3   attN;
    GA_RWHandleV3   attV;
    GA_RWHandleV3   attCenterUV;
    GA_RWHandleI    attId;
    GA_RWHandleF    attLife;
    GA_RWHandleI    attSpawn;
    GA_RWHandleI    attActive;
    GA_RWHandleI    attIsMature;
    GA_RWHandleI    attDensity;
    GA_RWHandleF    attBlend;
    GA_RWHandleF    attRandT;
    GA_RWHandleF    attMaxDeltaOnD;
    GA_RWHandleI    attDeleteFaster;
    GA_RWHandleV3   refAttV;
    GA_RWHandleV3   refAttN;
    GA_RWHandleV3   tangeant;
    GA_RWHandleI    attFadeIn;
    GA_RWHandleF    temporalComponentKt;
    GA_RWHandleV3   AttCd;
    GA_RWHandleI    attNumberOfPrimitives;
    GA_RWHandleI    attStatus;

    GA_RWHandleV3 attVSurface;
    //GA_RWHandleV3 attNSurface;
    GA_RWHandleF attDivergence;


};

}

#endif
