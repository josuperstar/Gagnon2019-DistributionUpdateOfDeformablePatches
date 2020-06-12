#include "HoudiniUtils.h"

#include <GU/GU_RandomPoint.h>


//============================== GET POINT NEIGHBORHOODS ======================================
using namespace TexturingFluids;
using namespace std;

std::set<GA_Offset> HoudiniUtils::GetNeighbors(GU_Detail *gdp,GA_Offset ppt)
{
    std::set<GA_Offset> results;

    GA_OffsetArray primitives;
    GA_Offset prim_off;

    GA_Size size = gdp->getPrimitivesReferencingPoint(primitives,ppt);


    for(GA_OffsetArray::const_iterator prims_it = primitives.begin(); prims_it != primitives.end(); ++prims_it)
    {
        prim_off = *prims_it;
        GA_Primitive *prim = gdp->getPrimitive(prim_off);
        //GEO_Primitive *geoPrim = gdp->getGEOPrimitive(prim_off);

        GA_Size numVertex = prim->getVertexCount();
        for (GA_Size i = 0; i< numVertex; i++)
        {
            GA_Offset point = prim->getPointOffset(i);
            if (point != ppt)
            {
                results.insert(point);
            }
        }
    }
    return results;

}

std::set<GA_Offset> HoudiniUtils::GetPrimitivesNeighbors(GU_Detail *gdp, GA_Offset ppt)
{
    std::set<GA_Offset> results;

    GA_OffsetArray primitives;
    GA_Offset prim_off;

    GA_Size size = gdp->getPrimitivesReferencingPoint(primitives,ppt);

    for(GA_OffsetArray::const_iterator prims_it = primitives.begin(); prims_it != primitives.end(); ++prims_it)
    {
        prim_off = *prims_it;
        results.insert(prim_off);
    }
    return results;
}


void HoudiniUtils::AttributeTransfert(GU_Detail *gdp, GU_Detail *referenceGdp, GEO_PointTreeGAOffset &tree, GA_PointGroup *grp, float transfertRadius, GA_RWHandleV3 refAtt, GA_RWHandleV3 destAtt)
{
    //============================== ATTRIB TRANSFERT ============================



    UT_Vector3 p;
    UT_Vector3 pN;

    UT_Vector3 nV; //neighbor value

    float W_k = 0.0f;
    float sumW_k = 0.0f;
    float distance = 0;
    //GA_PointGroup *grp = (GA_PointGroup *)gdp->pointGroups().find(this->markerGroupName.c_str());

    UT_Vector3 value;
    GA_Offset ppt;
    {
        // Transfert velocity from
        GA_FOR_ALL_GROUP_PTOFF(gdp,grp,ppt)
        {
            //get adjacent marker
            p = gdp->getPos3(ppt);
            value = UT_Vector3(0,0,0);

            GEO_PointTreeGAOffset::IdxArrayType close_particles_indices;

            tree.findAllCloseIdx(p,
                                 transfertRadius,
                                 close_particles_indices);

            unsigned close_particles_count = close_particles_indices.entries();
            if (close_particles_count > 0)
            {
                bool closeEnough = false;
                for(int j=0; j<close_particles_count;j++ )
                {
                    int index = close_particles_indices.array()[j];
                    pN = referenceGdp->getPos3(index);
                    nV = refAtt.get(index);

                    distance = distance3d(p,pN);
                    if (distance == 0)
                    {
                        W_k = 1;
                        closeEnough = true;
                    }
                    else
                    {
                        W_k = 1.f / (distance);//*distance*distance);
                    }
                    value += nV * W_k;
                    sumW_k += W_k;
                    if (closeEnough)
                        break;
                }
            }
            if (sumW_k <= 0)
            {
                //cout << "No neighborhood ?"<<endl;
                value = UT_Vector3(0,0,0);
            }
            else
            {
                value /= sumW_k;
                destAtt.set(ppt,value);
            }
            sumW_k = 0;
        }
    }

    //===============================================================================
}


void HoudiniUtils::AttributeTransfert(GU_Detail *gdp, GU_Detail *referenceGdp, GEO_PointTreeGAOffset &tree, GA_PointGroup *grp, float transfertRadius, GA_RWHandleV3 refAtt, GA_RWHandleV3 destAtt, float excludingDot)
{
    //============================== ATTRIB TRANSFERT ============================

    UT_Vector3 p;
    UT_Vector3 pN;


    UT_Vector3 nV; //neighbor value

    float W_k = 0.0f;
    float sumW_k = 0.0f;
    float distance = 0;
    UT_Vector3 N;
    //GA_PointGroup *grp = (GA_PointGroup *)gdp->pointGroups().find(this->markerGroupName.c_str());

    GA_RWHandleV3 attN(gdp->findFloatTuple(GA_ATTRIB_POINT,"N", 3));

    UT_Vector3 value;
    GA_Offset ppt;
    {
        // Transfert velocity from
        GA_FOR_ALL_GROUP_PTOFF(gdp,grp,ppt)
        {
            //get adjacent marker
            p = gdp->getPos3(ppt);
            value = UT_Vector3(0,0,0);
            N = attN.get(ppt);

            GEO_PointTreeGAOffset::IdxArrayType close_particles_indices;

            tree.findAllCloseIdx(p,
                                 transfertRadius,
                                 close_particles_indices);

            unsigned close_particles_count = close_particles_indices.entries();
            if (close_particles_count > 0)
            {
                bool closeEnough = false;
                for(int j=0; j<close_particles_count;j++ )
                {
                    int index = close_particles_indices.array()[j];
                    pN = referenceGdp->getPos3(index);
                    nV = refAtt.get(index);

                    if (dot(nV,N) < excludingDot)
                        continue;


                    distance = distance3d(p,pN);
                    if (distance == 0)
                    {
                        W_k = 1;
                        closeEnough = true;
                    }
                    else
                    {
                        W_k = 1.f / (distance);//*distance*distance);
                    }
                    value += nV * W_k;
                    sumW_k += W_k;
                    if (closeEnough)
                        break;
                }
            }
            if (sumW_k <= 0)
            {
                //cout << "No neighborhood ?"<<endl;
                value = UT_Vector3(0,0,0);
            }
            else
            {
                value /= sumW_k;
                destAtt.set(ppt,value);
            }
            sumW_k = 0;
        }
    }

    //===============================================================================
}


void HoudiniUtils::AttributeTransfert(GU_Detail *gdp, GEO_PointTreeGAOffset &tree, GA_PointGroup *fromGroup, GA_PointGroup *destGroup, float transfertRadius,
                                      GA_RWHandleV3 refAtt,
                                      GA_Attribute        *uvsAtt,
                                      const GA_AIFNumericArray *uvsArray)
{
    UT_Vector3 p;
    UT_Vector3 pN;

    UT_Vector3 nV; //neighbor value

    UT_FloatArray         uvArrayData;

    float W_k = 0.0f;
    float sumW_k = 0.0f;
    float distance = 0;
    //GA_PointGroup *grp = (GA_PointGroup *)gdp->pointGroups().find(this->markerGroupName.c_str());

    UT_Vector3 value;
    GA_Offset ppt;
    {
        // Transfert velocity from
        GA_FOR_ALL_GROUP_PTOFF(gdp,destGroup,ppt)
        {
            //get adjacent marker
            p = gdp->getPos3(ppt);
            value = UT_Vector3(0,0,0);

            GEO_PointTreeGAOffset::IdxArrayType close_particles_indices;

            tree.findAllCloseIdx(p,
                                 transfertRadius,
                                 close_particles_indices);

            unsigned close_particles_count = close_particles_indices.entries();
            if (close_particles_count > 0)
            {
                bool closeEnough = false;
                for(int j=0; j<close_particles_count;j++ )
                {
                    int index = close_particles_indices.array()[j];
                    if (!fromGroup->containsOffset(index))
                        continue;

                    pN = gdp->getPos3(index);
                    nV = refAtt.get(index);

                    distance = distance3d(p,pN);
                    if (distance == 0)
                    {
                        W_k = 1;
                        closeEnough = true;
                    }
                    else
                    {
                        W_k = 1.f / (distance);//*distance*distance);
                    }
                    value += nV * W_k;
                    sumW_k += W_k;
                    if (closeEnough)
                        break;
                }
            }
            if (sumW_k <= 0)
            {
                //cout << "No neighborhood ?"<<endl;
                value = UT_Vector3(0,0,0);
            }
            else
            {
                value /= sumW_k;
                //destAtt.set(ppt,value);
            }


            uvsArray->get(uvsAtt, ppt, uvArrayData);
            uvArrayData.append(value.x());
            uvArrayData.append(value.y());
            uvArrayData.append(value.z());
            // Write back
            uvsArray->set(uvsAtt, ppt, uvArrayData);



            sumW_k = 0;
        }
    }
}


void HoudiniUtils::AttributeTransfert(GU_Detail *gdp, GEO_PointTreeGAOffset &tree, GA_PointGroup *fromGroup, GA_PointGroup *destGroup, float transfertRadius,
                                      GA_RWHandleF refAtt,
                                      GA_Attribute        *uvsAtt,
                                      const GA_AIFNumericArray *uvsArray)
{
    UT_Vector3 p;
    UT_Vector3 pN;

    float nV; //neighbor value
    UT_FloatArray         uvArrayData;
    float W_k = 0.0f;
    float sumW_k = 0.0f;
    float distance = 0;
    //GA_PointGroup *grp = (GA_PointGroup *)gdp->pointGroups().find(this->markerGroupName.c_str());

    float value;
    GA_Offset ppt;
    {
        // Transfert velocity from
        GA_FOR_ALL_GROUP_PTOFF(gdp,destGroup,ppt)
        {
            //get adjacent marker
            p = gdp->getPos3(ppt);
            value = 0;

            GEO_PointTreeGAOffset::IdxArrayType close_particles_indices;

            tree.findAllCloseIdx(p,
                                 transfertRadius,
                                 close_particles_indices);

            unsigned close_particles_count = close_particles_indices.entries();
            if (close_particles_count > 0)
            {
                bool closeEnough = false;
                for(int j=0; j<close_particles_count;j++ )
                {
                    int index = close_particles_indices.array()[j];
                    if (!fromGroup->containsOffset(index))
                        continue;

                    pN = gdp->getPos3(index);
                    nV = refAtt.get(index);

                    distance = distance3d(p,pN);
                    if (distance == 0)
                    {
                        W_k = 1;
                        closeEnough = true;
                    }
                    else
                    {
                        W_k = 1.f / (distance);//*distance*distance);
                    }
                    value += nV * W_k;
                    sumW_k += W_k;
                    if (closeEnough)
                        break;
                }
            }
            if (sumW_k <= 0)
            {
                //cout << "No neighborhood ?"<<endl;
                value = 0;
            }
            else
            {
                value /= sumW_k;


                uvsArray->get(uvsAtt, ppt, uvArrayData);
                uvArrayData.append(value);
                // Write back
                uvsArray->set(uvsAtt, ppt, uvArrayData);
            }


            sumW_k = 0;
        }
    }
}


void HoudiniUtils::AttributeTransfert(GU_Detail *gdp, GU_Detail *referenceGdp, GEO_PointTreeGAOffset &tree, float transfertRadius, GA_RWHandleI refAtt, GA_RWHandleI destAtt)
{
    UT_Vector3 p;
    UT_Vector3 pN;


    int nV; //neighbor value

    float distance = 0;
    UT_Vector3 N;
    //GA_PointGroup *grp = (GA_PointGroup *)gdp->pointGroups().find(this->markerGroupName.c_str());

    GA_RWHandleI attActive(gdp->findIntTuple(GA_ATTRIB_POINT,"active", 1));
    int active = 0;
    GA_RWHandleV3 attN(gdp->findFloatTuple(GA_ATTRIB_POINT,"N", 3));

    int value;
    GA_Offset ppt;
    {
        // Transfert velocity from
        GA_FOR_ALL_PTOFF(gdp,ppt)
        {

            active = attActive.get(ppt);
            if (active == 0)
                continue;
            //get adjacent marker
            p = gdp->getPos3(ppt);
            value = refAtt.get(ppt);
            N = attN.get(ppt);

            nV = refAtt.get(ppt);
            GEO_PointTreeGAOffset::IdxArrayType close_particles_indices;

            tree.findAllCloseIdx(p,
                                 transfertRadius,
                                 close_particles_indices);

            unsigned close_particles_count = close_particles_indices.entries();
            if (close_particles_count > 0)
            {
                bool closeEnough = false;
                for(int j=0; j<close_particles_count; j++)
                {
                    int index = close_particles_indices.array()[j];

                    pN = referenceGdp->getPos3(index);


                    //if (dot(nV,N) < excludingDot)
                    //    continue;


                    distance = distance3d(p,pN);
                    if (distance == 0)
                    {

                        closeEnough = true;
                    }
                    if (transfertRadius > distance)
                        destAtt.set(index,1);

                }
            }

        }
    }

    //===============================================================================
}

UT_Vector3 HoudiniUtils::GetBarycentricPosition(UT_Vector3 A,UT_Vector3 B, UT_Vector3 C, UT_Vector3 a, UT_Vector3 b, UT_Vector3 c, UT_Vector3 position)
{
    double a1 = ((B.y()-C.y())*(position.x()-C.x())+(C.x()-B.x())*(position.y()-C.y()));
    double a2 = ((B.y()-C.y())*(A.x()-C.x())+(C.x()-B.x())*(A.y()-C.y()));
    double b1 = ((C.y()-A.y())*(position.x()-C.x())+(A.x()-C.x())*(position.y()-C.y()));
    double b2 = ((B.y()-C.y())*(A.x()-C.x())+(C.x()-B.x())*(A.y()-C.y()));

    if (a2 != 0 && b2 != 0)
    {
    double alpha = a1/a2 ;
        double beta = b1/b2 ;
        double gamma = 1-alpha-beta;
        a *= alpha;
        b *= beta;
        c *= gamma;

        return a+b+c;
    }
    return UT_Vector3{0,0,0};
}

