#include "PatchedSurfaceGagnon2019.h"
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
#include <GEO/GEO_PrimClassifier.h>
#include <GEO/GEO_PointClassifier.h>
#include <GEO/GEO_PrimConnector.h>
#include <GU/GU_NeighbourList.h>
#include <GU/GU_RayIntersect.h>
#include <GU/GU_Flatten.h>

//#include <Strategies/StrategySurfaceTextureSynthesis.h>
#include <Core/Gagnon2019/Bridson2012PoissonDiskDistribution.h>
#include <Core/HoudiniUtils.h>


PatchedSurfaceGagnon2019::PatchedSurfaceGagnon2019(GU_Detail *surface, GU_Detail *trackersGdp) : DeformableGridsManagerGagnon2019(surface, trackersGdp)
{
    //this->numberOfPatches = 0;

    this->maxId = 0;
    this->gridCenterPosition.clear();

    /*
    uvsArray->clear(this->uvsAtt);
    patchIdsAtt->clear(patchIdsArrayAttrib);
    alphaAtt->clear(alphaArrayAtt);
    */
    //=========================== PATCH ID ARRAY ATTRIB ==========================
    UT_String patchName(this->patchIdsName);
    patchIdsArrayAttrib = surface->findIntArray(GA_ATTRIB_POINT,
                                            patchName,
                                            // Allow any tuple size to match
                                            -1, -1);
    if (!patchIdsArrayAttrib)
    {
        patchIdsArrayAttrib = surface->addIntArray(GA_ATTRIB_POINT,
                                    patchName,
                    // Tuple size one means we group the array into
                    // logical groups of 1.  It does *NOT* affect
                    // the length of the arrays, which are always
                    // measured in ints.
                                    1);
    }
    // Make sure we are an array.  Note tuples do not match to this,
    // nor do GA_*Handle* match!
    // We will match both int and float here, however.
    // (For string, getAIFSharedStringArray)
    patchIdsAtt = patchIdsArrayAttrib->getAIFNumericArray();
    if (!patchIdsAtt)
    {
        UT_WorkBuffer           buf;
        buf.sprintf("Attribute \"%s\" not a numeric array!",
                    (const char *) patchName);
        //addError(SOP_MESSAGE, buf.buffer());
        //return error();
    }
    patchIdsArrayAttrib->clearDataId();

    //=========================== ALPHA ARRAY ATTRIB ==========================
    UT_String uvname(this->uvArrayName);
    uvsAtt = surface->findFloatArray(GA_ATTRIB_POINT,
                                            uvname,
                                            // Allow any tuple size to match
                                            -1, -1);

    if (!uvsAtt)
    {
        uvsAtt = surface->addFloatArray(GA_ATTRIB_POINT,
                                    uvname,
                    // Tuple size one means we group the array into
                    // logical groups of 1.  It does *NOT* affect
                    // the length of the arrays, which are always
                    // measured in ints.
                                    1);
    }

    uvsArray = uvsAtt->getAIFNumericArray();
    if (!uvsArray)
    {
        UT_WorkBuffer           buf;
        buf.sprintf("Attribute \"%s\" not a numeric array!",
                    (const char *) uvname);
        //addError(SOP_MESSAGE, buf.buffer());
        //return error();
    }
    uvsAtt->clearDataId();
    //=========================== ALPHA ARRAY ATTRIB ==========================

    UT_String aname(this->alphaArrayName);
    alphaArrayAtt = surface->findFloatArray(GA_ATTRIB_POINT,
                                            aname,
                                            // Allow any tuple size to match
                                            -1, -1);
    if (!alphaArrayAtt)
    {
        alphaArrayAtt = surface->addFloatArray(GA_ATTRIB_POINT,
                                    aname,
                    // Tuple size one means we group the array into
                    // logical groups of 1.  It does *NOT* affect
                    // the length of the arrays, which are always
                    // measured in ints.
                                    1);
    }
    // Make sure we are an array.  Note tuples do not match to this,
    // nor do GA_*Handle* match!
    // We will match both int and float here, however.
    // (For string, getAIFSharedStringArray)
    alphaAtt = alphaArrayAtt->getAIFNumericArray();
    if (!alphaAtt)
    {
        UT_WorkBuffer           buf;
        buf.sprintf("Attribute \"%s\" not a numeric array!",
                    (const char *) aname);
        //addError(SOP_MESSAGE, buf.buffer());
        //return error();
    }
    alphaArrayAtt->clearDataId();

    this->uvFlatteningTime = 0;
    this->gridMeshCreation = 0;
    this->gridAdvectionTime = 0;
    this->markerAdvectionTime = 0;
    this->patchCreationTime = 0;
    this->nbOfFlattenedPatch = 0;
    this->updatePatchesTime = 0;

}

PatchedSurfaceGagnon2019::~PatchedSurfaceGagnon2019()
{
    this->rays.clear();
}


//================================================================================================

//                                     POISSON DISK SAMPLING

//================================================================================================


void PatchedSurfaceGagnon2019::PoissonDiskSampling(GU_Detail *levelSet, GU_Detail *trackersGdp, ParametersDeformablePatches params)
{

    //This is a function that does a Poisson Disk Sampling using the approach of Bridson 2012 paper
    //This function is a wrapper to the Bridson2012PoissonDiskDistribution class.
    //It basically take the points from Houdini and fill it to the approach.


    std::clock_t addPoissonDisk;
    addPoissonDisk = std::clock();


    GEO_PointTreeGAOffset trackerTree;
    trackerTree.build(trackersGdp, NULL);
    int numberOfInitialTrackers = this->numberOfPatches;
    cout << this->approachName<<" initial number of point "<<numberOfInitialTrackers<< " "<<endl;
    //cout << "[Yu2011] we have "<<numberOfPoints << " existing point(s) in trackersGdp"<<endl;
    Bridson2012PoissonDiskDistributionGagnon2019 poissonDiskDistribution;
    poissonDiskDistribution.PoissonDiskSampling(trackersGdp, trackerTree, levelSet,params.poissondiskradius, params.poissonAngleNormalThreshold, params);
    this->numberOfNewPatches = poissonDiskDistribution.numberOfNewPoints;
    cout << this->approachName<<" poisson disk sample result: "<< this->numberOfNewPatches<< " new point(s)"<<endl;
    this->poissondisk += (std::clock() - addPoissonDisk) / (double) CLOCKS_PER_SEC;
    cout << this->approachName<<" Total :"<<trackersGdp->getNumPoints()<<endl;

}

GA_Offset PatchedSurfaceGagnon2019::CreateAPatch(GU_Detail *trackersGdp, UT_Vector3 position, UT_Vector3 normal, ParametersDeformablePatches params)
{
    //This is a function that does a Poisson Disk Sampling using the approach of Bridson 2012 paper
    //This function is a wrapper to the Bridson2012PoissonDiskDistribution class.
    //It basically take the points from Houdini and fill it to the approach.

    std::clock_t addPoissonDisk;
    addPoissonDisk = std::clock();


    GEO_PointTreeGAOffset trackerTree;
    trackerTree.build(trackersGdp, NULL);

    //cout << "[Yu2011] we have "<<numberOfPoints << " existing point(s) in trackersGdp"<<endl;
    Bridson2012PoissonDiskDistributionGagnon2019 poissonDiskDistribution;
    int numberOfClosePoint = 0;
    GA_Offset newPoint = poissonDiskDistribution.CreateAParticle(trackersGdp, trackerTree, position,normal,1,numberOfClosePoint,params);

    cout << this->approachName<<" poisson disk sample "<<trackersGdp->getNumPoints()<< " point(s)"<<endl;
    //this->numberOfPatches = trackersGdp->getNumPoints();

    this->poissondisk += (std::clock() - addPoissonDisk) / (double) CLOCKS_PER_SEC;
    return newPoint;
}



//================================================================================================

//                                       ADD PATCHES USING BARYCENTRIC COORDONATE

//================================================================================================


void PatchedSurfaceGagnon2019::AddDeformablePatchesUsingBarycentricCoordinates(GU_Detail *deformableGridsGdp,GU_Detail *surfaceGdp, GU_Detail *trackersGdp, ParametersDeformablePatches params, GEO_PointTreeGAOffset &surfaceTree,  GU_RayIntersect &ray)
{

    //This function is used to transfer the uv list from the deformable patches to the surface where the texture will be synthesis.
    cout << this->approachName<<"[AddPatchesUsingBarycentricCoordinates]" << endl;

    std::clock_t addPatchesStart;
    addPatchesStart = std::clock();


    float beta = params.Yu2011Beta;
    float d = params.poissondiskradius;
    float gridwidth = (2+beta)*d; //same formula used in DeformableGrids.cpp
    fpreal patchRadius = (fpreal)gridwidth;

    //================================ CREATE PATCH GROUPS ==============================
    //GA_PointGroup *grpMarker = (GA_PointGroup *)trackersGdp->pointGroups().find(this->markerGroupName.c_str());

    GA_GroupType primitiveGroupType = GA_GROUP_PRIMITIVE;
    const GA_GroupTable *primitiveGTable = deformableGridsGdp->getGroupTable(primitiveGroupType);

    /*
    GA_RWHandleV3 attN(trackersGdp->addFloatTuple(GA_ATTRIB_POINT,"N", 3));
    GA_RWHandleI attId(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"id",1));
    GA_RWHandleI attActive(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"active",1));
    GA_RWHandleF attLife(trackersGdp->findFloatTuple(GA_ATTRIB_POINT,"life",1));
    */
    GA_RWHandleV3 attNSurface(surfaceGdp->addFloatTuple(GA_ATTRIB_POINT,"N", 3));
    GA_RWHandleI attNumberOfPatch(surfaceGdp->addIntTuple(GA_ATTRIB_POINT,"numberOfPatch",1));

    if (attLife.isInvalid())
    {
        cout << this->approachName<<"particles have no life"<<endl;
        return;
    }


    set<int> patchTreated;
    float thresholdDistance = params.maximumProjectionDistance;
    float cs = params.CellSize;
    float padding = beta;
    float r = params.poissondiskradius*2+padding;
    float kd = params.poissondiskradius;
    GA_Offset ppt;
    UT_Vector3 N;
    UT_Vector3 NN;
    UT_Vector3 position;
    UT_Vector3 patchP;

    GA_Offset surfacePointOffset;
    int patchNumber = 0;
    {
        GA_FOR_ALL_PTOFF(trackersGdp, ppt)
        {
            patchNumber = attId.get(ppt);
            int active = attActive.get(ppt);
            float life = attLife.get(ppt);
            if (active == 0 && life <= 0 )
                continue;
            if (params.testPatch == 1 && params.patchNumber != patchNumber)
                continue;
            N = attN.get(ppt);
            position = trackersGdp->getPos3(ppt);

            UT_IntArray         patchArrayData;
            //getting neigborhood
            // Close particles indices
            GEO_PointTreeGAOffset::IdxArrayType surfaceNeighborhoodVertices;
            surfaceTree.findAllCloseIdx(position,
                                 patchRadius*2,
                                 surfaceNeighborhoodVertices);

            unsigned close_particles_count = surfaceNeighborhoodVertices.entries();

            string str = std::to_string(patchNumber);
            UT_String patchGroupName("patch"+str);
            //cout << "Create patch "<<patchGroupName<<endl;
            GA_PointGroup* patchGroup = surfaceGdp->newPointGroup(patchGroupName, 0);
            UT_String gridGroupName("grid"+str);
            GA_PrimitiveGroup*  gridPrimitiveGroup  = (GA_PrimitiveGroup*)primitiveGTable->find(gridGroupName);
            if (gridPrimitiveGroup == 0x0)
                continue;

            ray.init(deformableGridsGdp,gridPrimitiveGroup);

            set<GA_Offset> neighborhood;
            for(int j=0; j<close_particles_count;j++ )
            {
                surfacePointOffset = surfaceNeighborhoodVertices.array()[j];
                neighborhood.insert(surfacePointOffset);
            }

            set<GA_Offset>::iterator itG;
            for(itG = neighborhood.begin(); itG != neighborhood.end(); ++itG)
            {
                surfacePointOffset = *itG;
                NN = attNSurface.get(surfacePointOffset);
                float dotP = dot(N,NN); //exlude points that are not in the same plane.
                //if (dotP < params.angleNormalThreshold*2)
                //    continue;

                patchP = surfaceGdp->getPos3(surfacePointOffset);
                //respect poisson disk criterion
                //UT_Vector3 pos          = trackersGdp->getPos3(neighbor);
                UT_Vector3 pos          = patchP;
                //=====================================================
                // Put this section in a function, it also used in Poisson Disk Criteria

                bool outsideOfSmallEllipse = false;
                bool insideBigEllipse = false;

                Bridson2012PoissonDiskDistributionGagnon2019::IsInsideBlendingKernel(pos, position, N,r,cs, kd, outsideOfSmallEllipse, insideBigEllipse );

                if (!insideBigEllipse)
                    continue;
                //=====================================================

                // Surface vertex is inside the blending kernel
                // We now want to know if there is a grid primitive nearby

                //
                //We kill a particle when the distortion of its associated
                //grid becomes inconsistent with our quality criteria:
                //1) if the grid no longer covers the blending kernel,
                GU_MinInfo mininfo;
                mininfo.init(thresholdDistance*3,0.0001);
                ray.minimumPoint(patchP,mininfo);
                if (mininfo.prim == 0x0)
                {
                    continue;
                }

//                float u = mininfo.u1;
//                float v = mininfo.v1;

//                //get pos of hit
//                UT_Vector4 hitPos;
//                mininfo.prim->evaluateInteriorPoint(hitPos,mininfo.u1,mininfo.v1);
//                float dist = distance3d(patchP,hitPos);
//                if (dist > thresholdDistance)
//                {
//                    //we can't project the point on the deformable grid
//                    int active = attActive.get(ppt);
//                    if (active == 1)
//                    {
//                        cout << "the grid "<< patchNumber<<" no longer covers the blending kernel" <<endl;
//                        attActive.set(ppt,0);
//                    }
//                }

                patchGroup->addOffset(surfacePointOffset);
                patchIdsAtt->get(patchIdsArrayAttrib,surfacePointOffset, patchArrayData);
                int exist = patchArrayData.find(patchNumber);
                if (exist == -1)
                {
                    patchArrayData.append(patchNumber);

                    int numberOfPatch = attNumberOfPatch.get(surfacePointOffset);
                    numberOfPatch++;
                    attNumberOfPatch.set(surfacePointOffset,numberOfPatch);
                }
                patchIdsAtt->set(patchIdsArrayAttrib,surfacePointOffset, patchArrayData);
            }
            neighborhood.clear();
        }
    }

    {
        GA_FOR_ALL_PTOFF(trackersGdp, ppt)
        {
            patchNumber = attId.get(ppt);
            if (patchTreated.count(patchNumber) > 0)
            {
                cout << this->approachName<<"We already treated patch "<< patchNumber << endl;
                return;
            }
            patchTreated.insert(patchNumber);

            int active = attActive.get(ppt);
            float life = attLife.get(ppt);

            if (active == 0 && life <= 0 )
                continue;

            //for test purposes
            if (params.testPatch == 1 && params.patchNumber != patchNumber)
                continue;

            string str = std::to_string(patchNumber);
            UT_String pointGroupName("patch"+str);
            //cout << "Create primitive group" << str <<endl;
            GA_PrimitiveGroup* primGrp = surfaceGdp->newPrimitiveGroup(pointGroupName);
            GA_GroupType groupType = GA_GROUP_POINT;
            const GA_GroupTable *gtable = surfaceGdp->getGroupTable(groupType);
            GA_OffsetArray primitives;
            GA_PointGroup* pointGrp = (GA_PointGroup*)gtable->find(primGrp->getName());
            if (pointGrp == 0x0)
                continue;

            GA_FOR_ALL_GROUP_PTOFF(surfaceGdp,pointGrp,ppt)
            {
                surfaceGdp->getPrimitivesReferencingPoint(primitives,ppt);
                for(GA_OffsetArray::const_iterator prims_it = primitives.begin(); !prims_it.atEnd(); ++prims_it)
                {
                    primGrp->addOffset(*prims_it);
                }
            }
        }
     }
    this->patchCreationTime += (std::clock() - addPatchesStart) / (double) CLOCKS_PER_SEC;
}


//======================================================================================================================================
//                                                  UpdateUsingBridson2012PoissonDisk
//======================================================================================================================================

void PatchedSurfaceGagnon2019::DeleteUnusedPatches(GU_Detail *gdp, GU_Detail *trackersGdp, ParametersDeformablePatches params)
{
//    cout << this->approachName<<" Update Using Bridson 2012 Poisson Disk with "<<numberOfPatches<<" existing trackers"<<endl;
    std::clock_t startUpdatePatches;
    startUpdatePatches = std::clock();

    //--------------------------------------------------------------------------

    /*
    GA_RWHandleI attId(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"id",1));
    GA_RWHandleF attLife(trackersGdp->findFloatTuple(GA_ATTRIB_POINT,"life",1));
    GA_RWHandleI attSpawn(trackersGdp->findIntTuple(GA_ATTRIB_POINT,"spawn",1));
    GA_RWHandleI attActive(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"active", 1));

    GA_RWHandleI attIsMature(trackersGdp->findIntTuple(GA_ATTRIB_POINT,"isMature", 1));
    */
    GU_Detail::GA_DestroyPointMode mode = GU_Detail::GA_DESTROY_DEGENERATE;

    GA_Offset ppt;
    int beforeAddingNumber = this->numberOfPatches;
    this->numberOfPatches = 0;
    int deletedPatches = 0;

    GA_PointGroup *grpToDestroy = (GA_PointGroup *)trackersGdp->newPointGroup("ToDelete");
    GA_GroupType groupType = GA_GROUP_POINT;
    const GA_GroupTable *gtable = gdp->getGroupTable(groupType);
    GA_GroupType primGroupType = GA_GROUP_PRIMITIVE;
    const GA_GroupTable *gPrimTable = gdp->getGroupTable(primGroupType);
    set<int> toDelete;
    //--------------------------- DELETE DEAD PATCH --------------------------------------------
    {
        GA_FOR_ALL_PTOFF(trackersGdp,ppt)
        {
            int id = attId.get(ppt);
            float life = attLife.get(ppt);
            int active = attActive.get(ppt);
            if (active == 0 && life <= 0.0f && params.frame != params.startFrame)
            {
                //cout << "Deleting deformable grid "<<id<<" mature "<<attIsMature.get(ppt)<<endl;
                toDelete.insert(id);
                //numberOfPatches--;
                //numberOfConcealedPatches++;

                string str = std::to_string(id);
                string groupName = "grid"+str;
                GA_PrimitiveGroup *primGroup = (GA_PrimitiveGroup*)gPrimTable->find(groupName.c_str());
                GA_PointGroup* pointGrp = (GA_PointGroup*)gtable->find(groupName.c_str());

                if (pointGrp != 0x0)
                {
                    //cout << "delete points "<<endl;
                    gdp->deletePoints(*pointGrp,mode);
                    //cout << "delete point group "<<groupName<<endl;
                    gdp->destroyPointGroup(pointGrp);
                }
                if (primGroup != 0x0)
                {
                    //cout << "delete point group "<<groupName<<endl;
                    gdp->destroyPrimitiveGroup(primGroup);
                }
                deletedPatches++;
            }
            else
            {
                this->numberOfPatches++;
            }
        }
    }
    {
        GA_FOR_ALL_PTOFF(trackersGdp,ppt)
        {
            int id = attId.get(ppt);
            if (toDelete.count(id) > 0)
            {
                //cout << "Delete poisson disk point "<<id<<endl;
                grpToDestroy->addOffset(ppt);
            }
        }
    }
    //destroying trackers
    trackersGdp->deletePoints(*grpToDestroy,mode);
    trackersGdp->destroyPointGroup(grpToDestroy);

    //cout <<this->approachName<< " Added "<<(numberOfPatches-beforeAddingNumber) <<" new patches"<<endl;
    cout <<this->approachName<< " Removed "<<(deletedPatches)<<" dead patches "<<endl;
    this->updatePatchesTime += (std::clock() - startUpdatePatches) / (double) CLOCKS_PER_SEC;
    cout << this->approachName<<" TOTAL "<<numberOfPatches<< " patches"<<endl;
}

