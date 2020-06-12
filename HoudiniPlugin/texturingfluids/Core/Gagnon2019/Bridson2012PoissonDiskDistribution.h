#ifndef __Bridson2012PoissonDiskDistributionGagnon2019_h_
#define __Bridson2012PoissonDiskDistributionGagnon2019_h_
#include <vector>
#include <SOP/SOP_Node.h>
#include <GEO/GEO_PrimPart.h>
#include <GEO/GEO_PointTree.h>
//#include <GA_ElementGroup.h>
#include <vector>

#include <openvdb/openvdb.h>
#include <GU/GU_PrimVDB.h>
#include <openvdb/tools/Interpolation.h>
#include <Core/Deformations/ParametersDeformablePatches.h>

using namespace std;
namespace TexturingFluids {


//Based on paper from Bridson "Fast Poisson Disk Sampling in Arbritrary Dimensions"
// http://people.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf


/*
The algorithm takes as input the extent of the sample domain in
Rn, the minimum distance r between samples, and a constant k
as the limit of samples to choose before rejection in the algorithm
(typically k = 30).
*/

class Bridson2012PoissonDiskDistributionGagnon2019
{

public:

    //Bridson2012PoissonDiskDistribution(){}
    ~Bridson2012PoissonDiskDistributionGagnon2019()
    {
        cout << "[Bridson2012PoissonDiskDistribution] destrotying grid"<<endl;
        //backgroundGrid.~TreeDGrid();

    }
    void PoissonDiskSampling(GU_Detail* trackersGdp, GEO_PointTreeGAOffset &tree, GU_Detail *levelSet, float diskRadius, float angleNormalThreshold, ParametersDeformablePatches params);
    void SetNumberOfPoint(int data){this->numberOfPoints = data;}
    void initializeGrid(GEO_PointTreeGAOffset &tree, GU_Detail *trackerGdp, float diskRadius,  float angleNormalThreshold);
    void CreateAPointDisk(GU_Detail* trackersGdp, UT_Vector3 position, UT_Vector3 N);
    GA_Offset CreateAParticle(GU_Detail *trackerGdp, GEO_PointTreeGAOffset &tree, UT_Vector3 newPointPosition, UT_Vector3 newPointNormal, float killDistance, int &numberOfClosePoint, ParametersDeformablePatches &params);
    void SetMaxId(long data){maxId = data;}
    int numberOfNewPoints;

private:

    openvdb::Vec3f projectPointOnLevelSet(openvdb::Vec3f point, float distance, openvdb::Vec3f grad );

    bool RespectCriterion(GU_Detail* trackers, GEO_PointTreeGAOffset &tree, UT_Vector3 newPointPosition, UT_Vector3 newPointNormal,  float killDistance, int &numberOfClosePoint, GA_Offset exclude , ParametersDeformablePatches params);


    float poissonDiskRadius;    //radius
    int k;      //the limit of samples to choose before rejection in the algorithm, typically k = 30
    int numberOfPoints;


    int n = 3; // n-dimensional
    int t; // number of attemps
    float cellSize;
    bool Circle = true;

    long maxId = 0;


};
}

#endif
