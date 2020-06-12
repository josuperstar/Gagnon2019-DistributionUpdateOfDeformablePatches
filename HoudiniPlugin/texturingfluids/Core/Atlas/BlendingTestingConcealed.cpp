#include "BlendingTestingConcealed.h"
#include "../HoudiniUtils.h"

//================================= RASTERIZE PRIMITIVE =================================

Pixel BlendingTestingConcealed::Blend(GU_Detail* deformableGrids, int i, int j, float w, float h,
                                int pixelPositionX, int pixelPositionY,
                                vector<int> &sortedPatches,
                                vector<UT_Vector3> &surfaceUv,
                                vector<UT_Vector3> &surfacePosition,
                                map<int,UT_Vector3> &trackersNormal,
                                map<int,UT_Vector3> &trackersPosition,
                                map<int,UT_Vector3> &trackersUVPosition,
                                map<int, bool> &usePatches,
                                map<string,GU_RayIntersect*> &rays,
                                GA_RWHandleV3 &attPointUV,
                                map<int,float> &fading,
                                vector<ImageCV*> textureExemplars,

                                ParametersDeformablePatches params)
{

    //GA_RWHandleF attAlpha(deformableGrids->findFloatTuple(GA_ATTRIB_POINT,"Alpha",1));

    GA_GroupType primGroupType = GA_GROUP_PRIMITIVE;
    const GA_GroupTable *gPrimTable = deformableGrids->getGroupTable(primGroupType);

    float d = params.poissondiskradius;

    int tw = textureExemplars[0]->GetWidth();
    int th = textureExemplars[0]->GetHeight();

    float epsilon = 0.0001f;

    Pixel R_eq4 = Pixel(0,0,0);
    R_eq4.A = 1;
    Pixel color_wi_sum = Pixel(0,0,0);
    color_wi_sum.A = 1;


    Pixel color = Pixel(0,0,0);
    color.A = 1;

    Pixel Cf = Pixel(0,0,0);
    Cf.A = 1;

    float sumAlpha = 0.0f;


    //Equation 2, Quality of a triangle
    GA_RWHandleF    attQt(deformableGrids->findFloatTuple(GA_ATTRIB_PRIMITIVE,"Qt",1));
    GA_RWHandleF    attQv(deformableGrids->findFloatTuple(GA_ATTRIB_POINT,"Qv",1));
    GA_RWHandleI    attBorder(deformableGrids->findIntTuple(GA_ATTRIB_POINT,"border",1));

    UT_Vector3 pixelPositionOnSurface;

    //We don't work with an image with no width of height
    if (w-1.0f == 0.0f || h-1.0f == 0.0f)
        return color;

    pixelPositionOnSurface.x() = ((float)pixelPositionX/(w-1));
    pixelPositionOnSurface.y() = ((float)pixelPositionY/(h-1));
    pixelPositionOnSurface.z() = 0;

    UT_Vector3 positionOnSurface = HoudiniUtils::GetBarycentricPosition(surfaceUv[0],surfaceUv[1],surfaceUv[2],surfacePosition[0],surfacePosition[1],surfacePosition[2],pixelPositionOnSurface);

    float thresholdProjectionDistance = d/2.0f;

    //================================= SUMARIZE ========================================
    // compute color according to the list of patch
    //for each pixel patch loop

    //give more weight to the first patch and decrease afterward

    int k = sortedPatches.size();
    int numberOfHit = 0;
    int nbPatches = k;
    vector<int>::iterator itPatch;
    for(itPatch = --sortedPatches.end(); itPatch != --sortedPatches.begin(); itPatch--)
    {
        int patchId = *itPatch;

        if ( fading.find(patchId) == fading.end())
            continue;

        //-------------------------------------------------------------
        //get deformable grids according to patchId in deformableGrids

        string str = std::to_string(patchId);
        string groupName;
        groupName = "grid"+str;

        GU_RayIntersect *ray;

        map<string,GU_RayIntersect*>::const_iterator it = rays.find(groupName);
        if (it==rays.end())
            continue;
        ray = rays[groupName];

        //--------------------------------------------------
        //Can we project the pixel on the current patch ?
        //We may want to put this test outside this function, especially if we want to create a shader.
        GU_MinInfo mininfo;
        mininfo.init(thresholdProjectionDistance,0.0001f);

        ray->minimumPoint(positionOnSurface,mininfo);
        if (!mininfo.prim)
            continue;
        //--------------------------------------------------

        const GEO_Primitive *prim = mininfo.prim;

        UT_Vector3 primN = prim->computeNormal();
//        if (dot(primN, trackersNormal[patchId]) < 0.7)
//            continue;

        //We don't work with primitive that don't have at least 3 vertices
        if (prim->getVertexCount() < 3)
            continue;

        //------------------------------ PARAMETRIC COORDINATE -----------------------------------
        float u = mininfo.u1;
        float v = mininfo.v1;

        //get pos of hit
        UT_Vector4 hitPos;
        mininfo.prim->evaluateInteriorPoint(hitPos,mininfo.u1,mininfo.v1);
        float dist = distance3d(positionOnSurface,hitPos);
        if (dist > thresholdProjectionDistance)
            continue;

        if (dist > d/10.0f)

        {
            UT_Vector3 AB = (positionOnSurface - hitPos);
            AB = AB.normalize();

            float angle = dot(AB,primN);
            if (abs(angle) < 0.8)
                continue;
        }

        // We want to check that the pixel is not outside of the uv patch radius transferred on the surface.
        UT_Vector3 tracker = trackersPosition[patchId];
        float gridwidth = (2+params.Yu2011Beta)*params.poissondiskradius;
        float d_P = distance3d(hitPos,tracker);
        if (d_P > gridwidth)
            continue;

        //cout << "Hit pos "<<hitPos<<endl;

        GA_Offset vertexOffset0 = prim->getVertexOffset(0);
        GA_Offset vertexOffset1 = prim->getVertexOffset(1);
        GA_Offset vertexOffset2 = prim->getVertexOffset(2);

        //------------------- secion 3.3.3 Estimating Grid Distortion --------------

        //------- position in polygon of the deformable grid --------
        if (attPointUV.isInvalid())
            continue;

        GA_Offset pointOffset0  = deformableGrids->vertexPoint(vertexOffset0);
        GA_Offset pointOffset1  = deformableGrids->vertexPoint(vertexOffset1);
        GA_Offset pointOffset2  = deformableGrids->vertexPoint(vertexOffset2);

        UT_Vector3 v0 = attPointUV.get(pointOffset0);
        UT_Vector3 v1 = attPointUV.get(pointOffset1);
        UT_Vector3 v2 = attPointUV.get(pointOffset2);


        //--------------------------Equation 6--------------------
        //temporal componet

        //from the Rapport de Chercher: https://hal.inria.fr/inria-00355827v4/document
        //Also available in the doc directory of this project: doc/yu2009SPTA.pdf

        //Temporal weights Finally we fade particles in and out at their creation and destruction,
        //using two weights Fin(t) and Fout(t). The fading period τ can be long (in our
        //implementation we used τ = 5 seconds). Fout is mainly used to force the fading out
        //of grids whose distortion would stop increasing. And if τ is longer than the particles
        //lifetime, Fin and Fout rule the weights of all particles and are thus renormalized at the end


        UT_Vector3 centerUV = trackersUVPosition[patchId];//UT_Vector3(0.5,0.5,0.0);

        float s = params.UVScaling;

        v0 = UT_Vector3(v0.x()-centerUV.x(),v0.y()-centerUV.y(),v0.z()-centerUV.z());
        v1 = UT_Vector3(v1.x()-centerUV.x(),v1.y()-centerUV.y(),v1.z()-centerUV.z());
        v2 = UT_Vector3(v2.x()-centerUV.x(),v2.y()-centerUV.y(),v2.z()-centerUV.z());

        v0 *= s;
        v1 *= s;
        v2 *= s;

        v0 = UT_Vector3(v0.x()+centerUV.x(),v0.y()+centerUV.y(),v0.z()+centerUV.z());
        v1 = UT_Vector3(v1.x()+centerUV.x(),v1.y()+centerUV.y(),v1.z()+centerUV.z());
        v2 = UT_Vector3(v2.x()+centerUV.x(),v2.y()+centerUV.y(),v2.z()+centerUV.z());

        UT_Vector3 positionInPolygon = v0+u*(v1-v0)+v*(v2-v0);

        //-----------------------------------
        //Q_v quality of the vertex, value from 0 to 1
        //The weights are computed for each vertex. During reconstruction, weights at arbitrary locations are interpolated
        //from vertices values.
        int   d_V1 = 1-attBorder.get(pointOffset0);
        int   d_V2 = 1-attBorder.get(pointOffset1);
        int   d_V3 = 1-attBorder.get(pointOffset2);

        float Q_t1 = attQv.get(pointOffset0);
        float Q_t2 = attQv.get(pointOffset1);
        float Q_t3 = attQv.get(pointOffset2);

        float Q_t = Q_t1+u*(Q_t2-Q_t1)+v*(Q_t3-Q_t1);
        float d_V = d_V1+u*(d_V2-d_V1)+v*(d_V3-d_V1);

        if (Q_t < 0.001)
            continue;

        float   Q_V = Q_t*d_V;

        // ------------------------ EUCLEDIAN DISTANCE ------------------------

        float minD = d;
        float maxD = gridwidth; //edge region

        //d_V =0 if V ∈ grid boundary 1 otherwise

        if (d_P > maxD)
            d_V = 0.0f;

        float C_s = 0.0f;
        if (d_P > minD && d_P <= maxD)
            C_s = 1-((d_P-minD)/(maxD-minD));
        if (d_P <= minD)
            C_s = 1.0f;

        //----------------------------------------------------------------

        float K_s = C_s*d_V*Q_V;

        //K_s should be between 0 and 1
        if (K_s < 0)
            K_s = 0;
        else if (K_s > 1.0f)
            K_s = 1.0f;

        //--------------------------Equation 5--------------------
        // section 3.4.1 Vertex Weights
        //The weight for each vertex is defined as the product of a spatial component and a temporal component

        if (K_s < epsilon)
            continue;

        //-----------------------------------------------------------------
        //getting the color from the texture exemplar
        int i2 = static_cast<int>(floor(positionInPolygon.x()*tw));
        int j2 = ((int)th-1)-static_cast<int>(floor((positionInPolygon.y())*th));
        int seamCarvingIndex = ((1-K_s) * params.NumberOfTextureSampleFrame);

        textureExemplars[seamCarvingIndex]->GetColor(i2,j2,0,color);
        sumAlpha += color.A;
        if (color.A == 0)
            continue;

        //cout << "Using this one !!"<<endl;
        usePatches[patchId] = true;
        numberOfHit++;
        if (K_s == 1)
            break;

    }
    //========================= END SUM ==================================

    if (numberOfHit == 0)
        Cf.A == 0;
    Cf.A = sumAlpha;

    return Cf;
}




