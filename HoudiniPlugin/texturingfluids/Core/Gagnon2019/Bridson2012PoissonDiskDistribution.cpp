#include "Bridson2012PoissonDiskDistribution.h"
//#include "HoudiniUtils.h"

#include <GU/GU_RandomPoint.h>
#include <GEO/GEO_PrimVDB.h>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/math/Transform.h>
#include <Core/Deformations/ParametersDeformablePatches.h>
using namespace TexturingFluids;
using namespace std;


//================================================================================================

//                                      POISSON DISK SAMPLING

//================================================================================================

/*
We assume the surface geometry is given as a signed distance function:
this permits fast projection of points to the surface. Pseudocode
is provided in Algorithm 1. In the outer loop we search
for “seed” sample points on the surface, checking every grid cell
that intersects the surface (i.e. where the level set changes sign) so
we don’t miss any components: in a cell we take up to t attempts,
projecting random points from the cell to the surface and stopping
when one satisfies the Poisson disk criterion, i.e. is at least distance
r from existing samples. Once we have a seed sample, we continue
sampling from it, taking a step of size e · r from the previous sample
along a random tangential direction d, again projecting to the
surface and checking the Poisson disk criterion. Parameters t = 30
and e = 1.085 worked well, but could be further tuned.
*/

//Input: Level set φ, radius r, # attempts t, extension e
//Output: Sample set S

/*
1: for all grid cells C where φ changes sign do
2:  for t attempts do
3:      Generate random point p in C
4:      Project p to surface of φ
5:      if p meets the Poisson Disk criterion in S then
6:          S ← S ∪ {p}
7:          Break
8:  if no point was found in C then
9:      Continue
10: while new samples are found do
11:     Generate random tangential direction d to surface at p
12:     q ← p + d · e · r
13:     Project q to surface of φ
14:     if q meets the Poisson Disk criterion in S then
15:         S ← S ∪ {q}
16:         p ← q
*/


void Bridson2012PoissonDiskDistributionGagnon2019::PoissonDiskSampling(GU_Detail* trackersGdp, GEO_PointTreeGAOffset &tree, GU_Detail *levelSet, float diskRadius, float angleNormalThreshold, ParametersDeformablePatches params)
{
    cout << "[Bridson2012PoissonDiskDistributionGagnon2019] on level set using a threshold of "<<angleNormalThreshold<<endl;

    GA_RWHandleV3   attN(trackersGdp->findFloatTuple(GA_ATTRIB_POINT,"N", 3));
    GA_RWHandleI    attActive(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"active", 1));
    GA_RWHandleI    attId(trackersGdp->findIntTuple(GA_ATTRIB_POINT,"id",1));
    GA_RWHandleI    attDensity(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"density", 1));
    GA_RWHandleI    isTangeantTracker(trackersGdp->findIntTuple(GA_ATTRIB_POINT,"isTrangeantTracker",1));

    this->numberOfNewPoints = 0;

    // Find first vdb primitive of input 0
    GEO_Primitive* prim;
    GEO_PrimVDB* phi = 0x0;
    GA_FOR_ALL_PRIMITIVES(levelSet, prim)
    {
        if (phi = dynamic_cast<GEO_PrimVDB*>(prim))
            break;
    }

    if (!phi || !phi->hasGrid())
    {
        cout << "[Bridson2012PoissonDiskDistribution] Input geometry 0 has no VDB grid!"<<endl;
        return;
    }

//    cout << "Grid name: " << phi->getGridName() << std::endl;
//    cout << "Storage type: " << phi->getStorageType() << ", " << phi->getTupleSize() << std::endl;
//    cout << "JSON: " << phi->getJSON() << std::endl;

    float a = 0.25; //promote this variable to the user interface
    this->poissonDiskRadius = diskRadius;
    float killDistance = (1-a)*diskRadius/2;

    cout << "[Bridson2012PoissonDiskDistribution] We have a valid vdb"<<endl;

    GA_RWHandleI    attDeleteFaster(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"deleteFaster", 1));
    GA_Offset ppt;
    // Only if we want to delete too close patches, which is not the case when we don't use fading in

    GA_FOR_ALL_PTOFF(trackersGdp,ppt)
    {
        int numberOfClosePoint;

        if (isTangeantTracker.isValid())
        {
            if (isTangeantTracker.get(ppt) == 1)
            {
                continue;
            }
        }
        UT_Vector3 pointPosition = trackersGdp->getPos3(ppt);
        UT_Vector3 pointNormal   = attN.get(ppt);
        if (attId.get(ppt) > this->maxId)
            this->maxId = attId.get(ppt);


        bool meetPoissonDiskCriterion = this->RespectCriterion(trackersGdp,tree, pointPosition, pointNormal, killDistance,  numberOfClosePoint, ppt, params);
        attDensity.set(ppt,numberOfClosePoint);

        if (attActive.get(ppt) == 0)
            continue;
        //If we have fading in, we are using 2019's approach
        //if (params.fadingIn == 1)
        {
            attActive.set(ppt,meetPoissonDiskCriterion);
        }
        if (!meetPoissonDiskCriterion)
        {
            //cout << "We should delete point "<<attId.get(ppt)<<", is in kill distance" <<killDistance<<endl;
        }
    }


    t = 30;

//    cout << "Grid name: " << phi->getGridName() << std::endl;
//    cout << "Storage type: " << phi->getStorageType() << ", " << phi->getTupleSize() << std::endl;
//    cout << "JSON: " << phi->getJSON() << std::endl;

    openvdb::GridBase::Ptr ptr = phi->getGridPtr();
    openvdb::FloatGrid::Ptr gridSurface = openvdb::gridPtrCast<openvdb::FloatGrid>(ptr);
    if(!gridSurface)
    {
        cout << "[Bridson2012PoissonDiskDistribution] Surface grid can't be converted in FloatGrid"<<endl;
        return;
    }
    if ((gridSurface->getGridClass() != openvdb::GRID_LEVEL_SET))
    {
        cout<< "[Bridson2012PoissonDiskDistribution] Surface grid is not a Level-set FloatGrid!"<<endl;
        return;
    }

    //=================================================================
    //                         OPEN VDB ACCESSORS
    //=================================================================

    // Create the gradient field
    openvdb::tools::Gradient<openvdb::FloatGrid> gradientOperator(*gridSurface);
    openvdb::VectorGrid::ConstPtr gridGradient = gradientOperator.process();
    //gridGradient->setName(ssGrad.str());
    openvdb::VectorGrid::ConstAccessor accessorGradient = gridGradient->getConstAccessor();
    openvdb::tools::GridSampler<openvdb::VectorGrid::ConstAccessor, openvdb::tools::BoxSampler>
            samplerGradient(accessorGradient, gridGradient->transform());
     // Get the sampler for the boundary grid (tri-linear filtering, in surface grid index space)
    openvdb::FloatGrid::Accessor accessorSurface = gridSurface->getAccessor();
    openvdb::tools::GridSampler<openvdb::FloatGrid::Accessor, openvdb::tools::BoxSampler>
            samplerSurface(accessorSurface, gridSurface->transform());

    //=================================================================
    // 1: for all grid cells C where φ changes sign do
    //=================================================================
    int nbOfCell = 0;

    //cout << "[Bridson2012PoissonDiskDistribution] Step 1: for all grid cells C where φ changes sign do"<<endl;
    for (openvdb::FloatGrid::ValueOnCIter gridCellIt = gridSurface->cbeginValueOn(); gridCellIt; ++gridCellIt)
    {
        float x = gridCellIt.getCoord().x();
        float y = gridCellIt.getCoord().y();
        float z = gridCellIt.getCoord().z();

        float offset = 0.5;

        openvdb::Vec3f cellPosition(x,y,z);
        openvdb::Vec3f worldCellPos = gridSurface->transform().indexToWorld(cellPosition);
        float boundaryDist = samplerSurface.wsSample(worldCellPos);
        //openvdb::Vec3f p    = it.getCoord();
        //if (boundaryDist <= 0.0)// && grad.length() > 0.0)
        {
            //if it is not close to the surface, continue
            if (abs(boundaryDist) > params.CellSize/2.0f) // We should use a threshold defined by the user
                continue;
            //=================================================================
            //2:  for t attempts do
            //=================================================================
            bool ableToInsertPoint = false;
            for(int i =0; i < t; i++)
            {
                //=================================================================
                //3:      Generate random point p in C
                //=================================================================
                int seed = i;
                //we want it to oscillate between -0.5 and 0.5

                srand(seed);
                float rx = (((double) rand()/(RAND_MAX)-offset));
                srand(seed+1);
                float ry = (((double) rand()/(RAND_MAX)-offset));
                srand(seed+2);
                float rz = (((double) rand()/(RAND_MAX)-offset));

                openvdb::Vec3f randomPosition(rx,ry,rz);
                randomPosition *= params.CellSize;

                openvdb::Vec3f p = worldCellPos+randomPosition;

                float newPointDistance = samplerSurface.wsSample(p);
                if (abs(newPointDistance) > params.poissondiskradius)
                {
                    //cout << "random point is outside of range"<<endl;
                    continue;
                }
                //cout << "abs(newPointDistance) > poissonDiskRadius"<<endl;
                //=================================================================
                //4:      Project p to surface of φ
                //=================================================================
                openvdb::Vec3f grad = samplerGradient.wsSample(p);
                if (grad.length() < 0.0001)
                    continue;
                openvdb::Vec3f poissonDisk = projectPointOnLevelSet(p,newPointDistance,grad);
                UT_Vector3 newPointPosition = UT_Vector3(poissonDisk.x(),poissonDisk.y(),poissonDisk.z());
                grad.normalize();
                UT_Vector3 newPointNormal = UT_Vector3(grad.x(),grad.y(),grad.z());

                //=================================================================
                //5:      if p meets the Poisson Disk criterion in S then
                //=================================================================
                int numberOfClosePoint;
                //cout << "Trying to fit "<<newPointPosition<<endl;
                bool meetPoissonDiskCriterion = this->RespectCriterion(trackersGdp,tree, newPointPosition, newPointNormal, poissonDiskRadius, numberOfClosePoint, -1, params);
                if (meetPoissonDiskCriterion)
                {
                    //=================================================================
                    //6:          S ← S ∪ {p}
                    //=================================================================
                    bool isValid = this->CreateAParticle(trackersGdp,tree, newPointPosition, newPointNormal, poissonDiskRadius , numberOfClosePoint, params);
                    if (isValid)
                    {
                        break;
                    }
                }
            }
            if (!ableToInsertPoint)
            {
                //cout << "after "<<t<<" attemps, there is no possible insertion."<<endl;
            }
        }
        nbOfCell++;
    }
    //cout << nbOfCell << " cells have been treated."<<endl;
    return;
}



//================================================================================================

//                                 PROJECT POINT ON LEVEL SET

//================================================================================================


openvdb::Vec3f Bridson2012PoissonDiskDistributionGagnon2019::projectPointOnLevelSet(openvdb::Vec3f point, float distance, openvdb::Vec3f grad)
{
    //get the norm of the gradient
    openvdb::Vec3f gradNorm = grad;
    gradNorm.normalize();

    //=================================================================
    //4:      Project p to surface of φ
    //=================================================================
    //projection
    //cout << "old p "<<p<<endl;
    //p = p - dist * (grad/gradNorm);
    point = point - distance * gradNorm;

    return point;
}


//================================================================================================

//                                      INSERT POISSON DISK

//================================================================================================


GA_Offset Bridson2012PoissonDiskDistributionGagnon2019::CreateAParticle(GU_Detail *trackersGdp, GEO_PointTreeGAOffset &tree, UT_Vector3 p, UT_Vector3 N,  float killDistance , int &numberOfClosePoint, ParametersDeformablePatches &params)
{

    GA_RWHandleV3   attN(trackersGdp->addFloatTuple(GA_ATTRIB_POINT,"N", 3));
    GA_RWHandleI    attActive(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"active", 1));
    GA_RWHandleI    attDensity(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"density", 1));
    GA_RWHandleI    attId(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"id",1));
    GA_RWHandleF    attLife(trackersGdp->addFloatTuple(GA_ATTRIB_POINT,"life",1));
    GA_RWHandleI    attSpawn(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"spawn",1));
    GA_RWHandleI    attIsMature(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"isMature", 1));
    GA_RWHandleF    attMaxDeltaOnD(trackersGdp->addFloatTuple(GA_ATTRIB_POINT,"maxDeltaOnD",1));
    GA_RWHandleF    attExistingLife(trackersGdp->addFloatTuple(GA_ATTRIB_POINT,"life", 1));
    GA_RWHandleI    attNumberOfPrimitives(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"numberOfPrimitives", 1));


    //---------- Tangeant Tracker -----------
    UT_Vector3 defaultDirection(1,0,0);
    UT_Vector3 S,T;
    float Tlenght = params.poissondiskradius/2.0f;
    GA_Offset tracker_offset;

    int divider = 1;
    if (params.useTangeantTracker == 1)
        divider = 2;
    if (trackersGdp->getNumPoints()/divider > this->maxId) //existing points
    {
        //cout << "New Max Id = "<<trackersGdp->getNumPoints()/divider<<endl;
        this->maxId = trackersGdp->getNumPoints()/divider;
    }
    int id = this->maxId+1;
    this->maxId = id;
    //cout << "New Max Id = "<<this->maxId<<endl;
    GA_Offset newPoint = trackersGdp->appendPoint();
    trackersGdp->setPos3(newPoint, p);
    attN.set(newPoint,N);
    attActive.set(newPoint,true);
    attDensity.set(newPoint,numberOfClosePoint);
    attId.set(newPoint,id);
    attSpawn.set(newPoint,0);
    attLife.set(newPoint,0.001f);
    attNumberOfPrimitives.set(newPoint,0);
    attIsMature.set(newPoint,0);
    attMaxDeltaOnD.set(newPoint,0);

    if (params.useTangeantTracker == 1)
    {
        GA_RWHandleI    isTangeantTracker(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"isTrangeantTracker",1));
        //---------- ADD TANGEANT TRACKER ----------
        cout << "Add Tangeant Tracker"<<endl;
        //put this in a function, and/or move this where we already add point
        S = cross(N,defaultDirection);
        S.normalize();
        T = cross(S,N);
        T.normalize();
        UT_Vector3 translation = T*Tlenght;
        //cout << "Translation: "<<translation<<endl;
        UT_Vector3 tangeantPosition = p + translation;
        //cout << "adding tangeant tracker"<<endl;
        tracker_offset = trackersGdp->appendPoint();
        isTangeantTracker.set(tracker_offset,1);
        attId.set(tracker_offset,id);
        trackersGdp->setPos3(tracker_offset,tangeantPosition);
    }

    if(params.startFrame == params.frame)
    {
        attExistingLife.set(newPoint,params.fadingTau);
    }
    tree.build(trackersGdp);

    this->numberOfNewPoints++;

    return newPoint;
}

void Bridson2012PoissonDiskDistributionGagnon2019::CreateAPointDisk(GU_Detail* trackersGdp, UT_Vector3 position, UT_Vector3 N)
{
    GA_RWHandleV3   attN(trackersGdp->addFloatTuple(GA_ATTRIB_POINT,"N", 3));
    GA_RWHandleI    attActive(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"active", 1));
    GA_RWHandleI    attDensity(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"density", 1));
    GA_RWHandleI    attId(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"id",1));
    GA_RWHandleF    attLife(trackersGdp->addFloatTuple(GA_ATTRIB_POINT,"life",1));
    GA_RWHandleI    attSpawn(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"spawn",1));
    GA_RWHandleI    attIsMature(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"isMature", 1));
    GA_RWHandleF    attMaxDeltaOnD(trackersGdp->addFloatTuple(GA_ATTRIB_POINT,"maxDeltaOnD",1));
    GA_RWHandleI    isTangeantTracker(trackersGdp->findIntTuple(GA_ATTRIB_POINT,"isTrangeantTracker",1));

    int divider = 1;
    if (isTangeantTracker.isValid())
        divider = 2;
    if (trackersGdp->getNumPoints()/divider > this->maxId) //existing points
    {
        cout << "New Max Id = "<<trackersGdp->getNumPoints()/divider<<endl;
        this->maxId = trackersGdp->getNumPoints()/divider;
    }
    int id = this->maxId+1;
    this->maxId = id;

    GA_Offset newPoint = trackersGdp->appendPoint();
    trackersGdp->setPos3(newPoint, position);
    attN.set(newPoint,N);
    attActive.set(newPoint,true);
    attDensity.set(newPoint,0);
    attId.set(newPoint,id);
    attSpawn.set(newPoint,0);
    attLife.set(newPoint,0.001f);
    attIsMature.set(newPoint,0);
    attMaxDeltaOnD.set(newPoint,0);

}

//================================================================================================

//                                      RESPECT CRITERION

//================================================================================================

bool Bridson2012PoissonDiskDistributionGagnon2019::RespectCriterion(GU_Detail* trackersGdp, GEO_PointTreeGAOffset &tree, UT_Vector3 newPointPosition, UT_Vector3 newPointNormal, float killDistance, int &numberOfClosePoint,   GA_Offset exclude, ParametersDeformablePatches params )
{
    numberOfClosePoint = 0;

    GA_RWHandleV3   attN(trackersGdp->addFloatTuple(GA_ATTRIB_POINT,"N", 3));
    GA_RWHandleI    attActive(trackersGdp->findIntTuple(GA_ATTRIB_POINT,"active", 1));

    GEO_PointTreeGAOffset::IdxArrayType close_particles_indices;

    tree.findAllCloseIdx(newPointPosition,
                         params.poissondiskradius*2,
                         close_particles_indices);

    int l = (int)close_particles_indices.entries();
    GA_Offset neighbor;
    bool tooClose = false;

    float cs    = params.CellSize;
    float r     = params.poissondiskradius;

    newPointNormal.normalize();
    float kd = killDistance;

    UT_Vector3 defaultDirection(1.012f,0.123f,0.002f);
    UT_Vector3 S,T;

    for(int j=0; j<l;j++)
    {
        neighbor = close_particles_indices.array()[j];
        if (attActive.get(neighbor) == 0)
            continue;
        if (neighbor == exclude)
            continue;

        UT_Vector3 pos          = trackersGdp->getPos3(neighbor);

        UT_Vector3 N            = attN.get(neighbor);
        N.normalize();
        S = cross(N,defaultDirection);
        S.normalize();
        T = cross(S,N);
        T.normalize();


        // Transform into local patch space (where STN is aligned with XYZ at the origin)
        const UT_Vector3 relativePosistion = pos - newPointPosition;
        UT_Vector3 poissonDiskSpace;
        poissonDiskSpace.x() = relativePosistion.dot(S);
        poissonDiskSpace.y() = relativePosistion.dot(T);
        poissonDiskSpace.z() = relativePosistion.dot(N);

        float dotN              = dot(N,newPointNormal);
        bool samePlane          = dotN > params.poissonAngleNormalThreshold;

        //(x/a)2 + (y/b)2 + (z/c)2 = 1
        float x = poissonDiskSpace.x();
        float y = poissonDiskSpace.y();
        float z = poissonDiskSpace.z();
        float a = r;
        float b = r;
        float c = cs*2;

        float a2 = kd;
        float b2 = kd;
        float c2 = cs;

        float smallEllipse = (x/a2)*(x/a2) + (y/b2)*(y/b2) + (z/c2)*(z/c2);
        float bigEllipse = (x/a)*(x/a) + (y/b)*(y/b) + (z/c)*(z/c);

        bool outsideOfSmallEllipse = false;
        bool insideBigEllipse = false;

        if (bigEllipse <= 1)
            insideBigEllipse = true;
        if (smallEllipse > 1)
            outsideOfSmallEllipse = true;

        //It is too close to the current point ?
        if(samePlane && !outsideOfSmallEllipse)
        {
            tooClose = true;
        }

        if(insideBigEllipse && samePlane)
            numberOfClosePoint++;
    }
    return !tooClose;
}
