#include "DistortionMetricSorkine2002.h"
#include <algorithm>    // std::max

bool DistortionMetricSorkine2002::ComputeDistortion(GU_Detail *trackersGdp, GU_Detail *deformableGridsGdp, GA_Offset trackerPpt,GA_PointGroup* pointGrp,GA_PrimitiveGroup *primGroup, ParametersDistortion params)
{
    //cout << "DISTORTION : Yu2011 Distortion using homogeneous removal blending"<<endl;
    //Sorkine 2002 distortion metric: https://igl.ethz.ch/projects/parameterization/BDPMP/sorkine02.pdf
    //Also available in the doc directory of this project

    GA_RWHandleV3    attSs(deformableGridsGdp->addFloatTuple(GA_ATTRIB_PRIMITIVE,"Ss",3));
    GA_RWHandleV3    attSt(deformableGridsGdp->addFloatTuple(GA_ATTRIB_PRIMITIVE,"St",3));
    GA_RWHandleF    attGMax(deformableGridsGdp->addFloatTuple(GA_ATTRIB_PRIMITIVE,"gmax",1));
    GA_RWHandleF    attGMin(deformableGridsGdp->addFloatTuple(GA_ATTRIB_PRIMITIVE,"gmin",1));
    GA_RWHandleF    attRefDistortion(deformableGridsGdp->findFloatTuple(GA_ATTRIB_PRIMITIVE,"distortion",1));
    GA_RWHandleF    attDt(deformableGridsGdp->addFloatTuple(GA_ATTRIB_PRIMITIVE,"dt",1));
    GA_RWHandleF    attQt(deformableGridsGdp->addFloatTuple(GA_ATTRIB_PRIMITIVE,"Qt",1));
    GA_RWHandleF    attQv(deformableGridsGdp->addFloatTuple(GA_ATTRIB_POINT,"Qv",1));

    GA_RWHandleV3   attRP(deformableGridsGdp->findFloatTuple(GA_ATTRIB_VERTEX,"refPosition",3));

    GA_RWHandleI    attId(trackersGdp->addIntTuple(GA_ATTRIB_POINT,"id",1));
    //GA_RWHandleF    attLife(trackersGdp->findFloatTuple(GA_ATTRIB_POINT,"life",1));
    GA_RWHandleI    attSpawn(trackersGdp->findIntTuple(GA_ATTRIB_POINT,"spawn",1));
    GA_RWHandleI    attActive(trackersGdp->findIntTuple(GA_ATTRIB_POINT,"active", 1));
    GA_RWHandleF    attMaxQt(trackersGdp->addFloatTuple(GA_ATTRIB_POINT,"maxQt", 1));



    bool killParticle = false;
    bool overDistorted = false;
    {
        GEO_Primitive *prim;
        float area;
        int nbVertex;

        float minQt = 10000.0f;
        GA_Offset primOffset;
        GA_FOR_ALL_GROUP_PRIMITIVES(deformableGridsGdp,primGroup,prim)
        {
            primOffset = prim->getMapOffset();
                    //move this into a function
            vector<UT_Vector3> trianglePoints;
            vector<GA_Offset> offsetArray;
            nbVertex = prim->getVertexCount();
            for (int i = 0; i < nbVertex; i++)
            {
                GA_Offset vertexPoint = prim->getVertexOffset(i);
                GA_Offset point = deformableGridsGdp->vertexPoint(vertexPoint);
                offsetArray.push_back(point);
                trianglePoints.push_back(deformableGridsGdp->getPos3(point));
            }

            GA_Offset vertexA = prim->getVertexOffset(0);
            GA_Offset vertexB = prim->getVertexOffset(1);
            GA_Offset vertexC = prim->getVertexOffset(2);

            UT_Vector3 p1 = attRP.get(vertexA);
            UT_Vector3 p2 = attRP.get(vertexB);
            UT_Vector3 p3 = attRP.get(vertexC);

            float s1 = p1.x();
            float s2 = p2.x();
            float s3 = p3.x();

            float t1 = p1.y();
            float t2 = p2.y();
            float t3 = p3.y();

            UT_Vector3 q1 = trianglePoints[0];
            UT_Vector3 q2 = trianglePoints[1];
            UT_Vector3 q3 = trianglePoints[2];

            //======================== SORKINE 2002 SECTION 3.2 =========================

            float dmax = params.Yu2011DMax;

            if (dmax < 1.0f)
                dmax = 1.1f;

            area = ((s2 - s1)*(t3-t1) - (s3-s1)*(t2-t1))/2;
            if (area != 0.0f)
            {
                UT_Vector3 Ss = ( (t2-t3)*q1 + (t3-t1) *q2+ (t1-t2)*q3)/(2*area);
                UT_Vector3 St = ( (s3-s2)*q1 + (s1-s3) *q2+ (s2-s1)*q3)/(2*area);

                float a = dot(Ss,Ss);
                float b = dot(Ss,St);
                float c = dot(St,St);

                float gmax = sqrt(0.5*((a+c)+sqrt((a-c)*(a-c) + 4*(b*b))));
                float gmin = sqrt(0.5*((a+c)-sqrt((a-c)*(a-c) + 4*(b*b))));

                //Yu 2011 Equation 1
                //Sorkine 2002: Note that dt >=1, and the equality holds if and only if T and T' are isometric
                float dt = std::max(gmax,1.0f/gmin);


                //Yu 2011 Equation 2
                //We then define the quality of a triangle as the ratio of its distortion with the maximum acceptable distortion,
                //δmax is then supposed to be greater than 1.
                //acording to the thesis https://tel.archives-ouvertes.fr/tel-00528781/document, page 74:
                //To limit the maximum distortion of an
                //advected texture, we kill a particle when D > δkill. We use δkill = 2.5 in our implementation.
                float Qt = std::max((dmax-dt)/(dmax-1.0f),0.0f);

                //Qt is equal to 1 for an un-distorted triangle, and is equal to 0 for a triangle where the distortion is larger
                //than δmax.

                //for debuging purpose:
                attSs.set(prim->getMapOffset(),Ss);
                attSt.set(prim->getMapOffset(),St);
                attGMax.set(prim->getMapOffset(),gmax);
                attGMin.set(prim->getMapOffset(),gmin);
                attDt.set(prim->getMapOffset(),dt);
                //The vertex quality measure. It is also use in the blending function.
                attQt.set(prim->getMapOffset(),Qt);

                //For each grid vertex V , we then compute its quality, QV as the mean of the quality of its incident triangles.
                //TODO compute Qv:

                //Apparently, this part does not work.
                //TODO: investigate

                GA_OffsetArray primitives;
                vector<GA_Offset>::iterator itPoint;
                GA_OffsetArray::const_iterator itPrim;
                float sumQt = 0;
                float nbQt = 0;
                for(itPoint = offsetArray.begin(); itPoint != offsetArray.end(); itPoint++)
                {
                    deformableGridsGdp->getPrimitivesReferencingPoint(primitives,*itPoint);
                    for(itPrim = primitives.begin(); itPrim != primitives.end(); itPrim++)
                    {
                        if (itPrim.item() != primOffset)
                        {
                            //this is an incident triangle:
                            float nqt = attQt.get(itPrim.item());
                            sumQt += nqt;
                            nbQt++;
                        }
                    }
                    float Qv = 0;
                    //cout << "sumQt "<<sumQt<<" nbQt "<<nbQt;
                    if (nbQt > 0)
                    {
                        Qv = sumQt/nbQt;
                        float currentQv = attQv.get(*itPoint);
                        if (currentQv > Qv)
                            attQv.set(*itPoint, Qv);
                    }
                    else
                    {
                        attQv.set(*itPoint, 0);
                    }
                }


                //We kill a particle if, for any vertex in the grid, we have QV < 1/2 (i.e., we keep a margin of quality for the fading-out).
                float QvMin = params.QvMin;
                if(Qt < QvMin)
                {
                    killParticle = true;
                }
                if (Qt <= 0.0f)
                {
                    overDistorted = true;
                }

                if (minQt > Qt)
                    minQt = Qt;

            }
            //====================================================================
            offsetArray.clear();
        }
        attMaxQt.set(trackerPpt,minQt);
    }
    if (killParticle)
    {
        //flag killing particle if it's not already flagged:
        if (params.flagDistortedParticles && attActive.get(trackerPpt) == 1)
        {
            attActive.set(trackerPpt,0);
            return true;
        }
    }
    return false;
}
