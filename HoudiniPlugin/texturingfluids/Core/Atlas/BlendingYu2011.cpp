#include "BlendingYu2011.h"
#include "../HoudiniUtils.h"

//================================= RASTERIZE PRIMITIVE =================================

Pixel BlendingYu2011::Blend(GU_Detail* deformableGrids, int i, int j, float w, float h,
                                int pixelPositionX, int pixelPositionY,
                                vector<int> &sortedPatches,
                                vector<UT_Vector3> &surfaceUv,
                                vector<UT_Vector3> &surfacePosition,
                                map<int,UT_Vector3> &trackersPosition,
                                map<int,UT_Vector3> &trackersUVPosition,
                                map<string,GU_RayIntersect*> &rays,
                                map<int,Pixel> &patchColors,
                                Pixel RM,           //Mean Value
                                GA_RWHandleV3 &attPointUV,
                                map<int,float> &fading,
                                vector<ImageCV*> textureExemplars,
                                ImageCV *displacementMapImage,
                                bool computeDisplacement,
                                bool renderColoredPatches,
                                Pixel &R_eq3,
                                Pixel &displacementSumEq3,
                                Pixel &displacementSumEq4,
                                ParametersDeformablePatches params)
{

    //GA_RWHandleF attAlpha(deformableGrids->findFloatTuple(GA_ATTRIB_POINT,"Alpha",1));
    bool useLocalRayIntersect = false;
    std::ofstream outfile;
    bool outputCSV = false;
    //if (outputCSV)
    //    cout << "----------------------"<<endl;

    GA_GroupType primGroupType = GA_GROUP_PRIMITIVE;
    const GA_GroupTable *gPrimTable = deformableGrids->getGroupTable(primGroupType);

    float d = params.poissondiskradius;
    Pixel displaceMean;
    if(computeDisplacement)
    {
        displaceMean = displacementMapImage->MeanValue();
    }

    int tw = textureExemplars[0]->GetWidth();
    int th = textureExemplars[0]->GetHeight();

    float epsilon = 0.0001f;

    Pixel R_eq4 = Pixel(0,0,0);
    R_eq4.A = 1;
    Pixel color_wi_sum = Pixel(0,0,0);
    color_wi_sum.A = 1;

    vector<Pixel> colorsList;
    vector<float> w_i_list;

    R_eq3 = Pixel(0,0,0);
    R_eq3.A = 1;
    float sumW2 = 0;
    float sumW = 0;

    Pixel color = Pixel(0,0,0);
    color.A = 1;

    Pixel displacement = Pixel(0,0,0);

    //Equation 2, Quality of a triangle
    GA_RWHandleF    attQt(deformableGrids->findFloatTuple(GA_ATTRIB_PRIMITIVE,"Qt",1));
    GA_RWHandleF    attQv(deformableGrids->findFloatTuple(GA_ATTRIB_POINT,"Qv",1));
    GA_RWHandleI    attBorder(deformableGrids->findIntTuple(GA_ATTRIB_POINT,"border",1));
    if (attBorder.isInvalid())
        return R_eq4;
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
        if (useLocalRayIntersect)
        {
            GA_PrimitiveGroup *primGroup = (GA_PrimitiveGroup*)gPrimTable->find(groupName.c_str());
            ray = new GU_RayIntersect(deformableGrids,primGroup);
            ray->init();
        }
        else
        {
            map<string,GU_RayIntersect*>::const_iterator it = rays.find(groupName);
            if (it==rays.end())
                continue;
            ray = rays[groupName];
        }
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
        //We don't work with primitive that don't have at least 3 vertices
        if (prim->getVertexCount() < 3)
            continue;

        //------------------------------ PARAMETRIC COORDINATE -----------------------------------
        float u = mininfo.u1;
        float v = mininfo.v1;

        //get pos of hit
        UT_Vector4 hitPos;
        mininfo.prim->evaluateInteriorPoint(hitPos,mininfo.u1,mininfo.v1);
        if (distance3d(positionOnSurface,hitPos) > thresholdProjectionDistance)
            continue;

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

        //The temporal component is simply a linear fade-in at the beginning of the life of a particle and a linear fade-out
        //after the particle has been killed.
        //Here, we take the fading from the particle stored in a map where the indexes are the patch number.
        float K_t = fading[patchId];
        if (K_t < 0.0f)
            K_t = 0.0f;
        else if (K_t > 1.0f)
            K_t = 1.0f;

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
        //float   Q_t = attQt.get(prim->getMapOffset());

        if (Q_t < 0.001)
            continue;

        float   Q_V = Q_t;

        //-----------------------------------------------------------------
        //getting the color from the texture exemplar
        int i2 = static_cast<int>(floor(positionInPolygon.x()*tw));
        int j2 = ((int)th-1)-static_cast<int>(floor((positionInPolygon.y())*th));

        //--------------------------Equation 7--------------------
        //spatial component

        /*
        The spatial component merges three factors: the quality around each grid vertex (QV , defined in section 3.3.3),
        a fall-off with the distance to the particle (in our implementation we take it linear), and a continuity factor
        ensuring a weight 0 on the boundary of the grid (to avoid spatial discontinuities during blending):

        K_s(V) = (1- (||v||-p)/d)*d_V*Q_V
        where d_V =0 if V ∈ grid boundary
        1 otherwise


        */

        //Here, we use a alpha chanel with linear fading from the center to compute the fall-off
        //Qv is an interpolation of the alpha chanel of the polygon use in the grid, for this patch.
        //Therefore, the Quality of the vertex has been computed before and stored in the alpha chanel

        //float dP = ((falloff.R+falloff.G+falloff.B)/3);


        float d_P = distance3d(positionInPolygon,centerUV);
        //float maxDUV = 0.175f; //should comme from the scaling used for the uv projection.
        //float maxDUV = (0.5f*sqrt(1.0f/params.UVScaling))/2.0f;

        float minDUV = 0.125*params.PatchScaling;
        float maxDUV = 0.25*params.PatchScaling; //blending region
        //float maxDUV = 0.5f;
        //d_V =0 if V ∈ grid boundary 1 otherwise
        //float d_V = 1.0f;

        //test
        //maxDUV = params.poissondiskradius/2;
        if (d_P > maxDUV)
            d_V = 0.0f;

        float C_s = 0.0f;
        if (d_P > minDUV && d_P <= maxDUV)
            C_s = 1-((d_P-minDUV)/(maxDUV-minDUV));
        if (d_P <= minDUV)
            C_s = 1.0f;

        //float K_s = (1.0f-(d_P/maxDUV))*d_V*Q_V;
        float K_s = C_s*d_V*Q_V;
        //cout << "Ks "<<K_s<<endl;
        //cout << "d_V "<<d_V<<endl;

        //K_s should be between 0 and 1
        if (K_s < 0)
            K_s = 0;
        else if (K_s > 1.0f)
            K_s = 1.0f;

        //--------------------------Equation 5--------------------
        // section 3.4.1 Vertex Weights
        //The weight for each vertex is defined as the product of a spatial component and a temporal component
        float w_v = K_s * K_t;
        if (w_v < epsilon)
            continue;



        //int seamCarvingIndex = ((w_v) * params.NumberOfTextureSampleFrame)-1;
        int textureSampleIndex = 0;
        if (renderColoredPatches)
            //set random colors per patch
            color = patchColors[patchId];
        else
        {
            textureExemplars[textureSampleIndex]->GetColor(i2,j2,0,color);
        }
            //cout << "Animated Color "<<color.R<<" "<<color.G<<" "<<color.B<<endl;

        if (computeDisplacement)
            displacementMapImage->GetColor(i2,j2,0,displacement);

        //clamping color values ...
        if (color.B > 1.0f)
            color.B = 1.0f;
        if (color.G > 1.0f)
            color.G = 1.0f;
        if (color.R > 1.0f)
            color.R = 1.0f;

        //cout << "w_i "<<w_v<<endl;
        w_i_list.push_back(w_v);

        //--------------------------------------------------------
        //Section 3.5.1
        //• Allocate the intermediate texture R at the required resolution (see discussion in 3.5.3), with one channel
        //for each aj (RGB), plus one channel for the wi (resp. wi2 if using Eq. 4).

        //• For each particle i: splat its grid into the texture (e.g., using a render target and the ordinary drawing API),
        //thus accumulating the sum(wi(x)aj (ui(x)) and the sum(wi(x)) (resp. sum(w2i) ) into their respective channels.

        sumW2 += w_v*w_v;
        sumW += w_v;



        colorsList.push_back(color);
        //cout <<color.R << "," <<color.G << "," <<color.B<<endl;
        //---------------- YU 2011 Equation 3 -----------------------
        //blending function
        R_eq3.R += w_v*(color.R);
        R_eq3.G += w_v*(color.G);
        R_eq3.B += w_v*(color.B);

        if (computeDisplacement)
        {
            displacementSumEq3.R += w_v*(displacement.R);
            displacementSumEq3.G += w_v*(displacement.G);
            displacementSumEq3.B += w_v*(displacement.B );
        }

        //---------------- YU 2011 Equation 4 -----------------------
        //blending function
        color_wi_sum.R += w_v*(color.R - RM.R);
        color_wi_sum.G += w_v*(color.G - RM.G);
        color_wi_sum.B += w_v*(color.B - RM.B);
        //blending function
        //R_eq4.R += (color.R - RM.R);
        //R_eq4.G += (color.G - RM.G);
        //R_eq4.B += (color.B - RM.B);


        //cout << "w_i "<< w_v<<endl;
        //cout << "color "<<color.R<< " "<<color.G<<" "<<color.B<<endl;

        if (computeDisplacement)
        {
            displacementSumEq4.R += w_v*(displacement.R - displaceMean.R);
            displacementSumEq4.G += w_v*(displacement.G - displaceMean.G);
            displacementSumEq4.B += w_v*(displacement.B - displaceMean.B);
        }
        k--;

        if (useLocalRayIntersect)
        {
            delete ray;
        }

    }
    //========================= END SUM ==================================


    //if the sumW is close to 0, that means we should have a new poisson disk there.
    //But how to handle the creation of a new patch from here ?

    if (sumW <= epsilon && sumW >= -epsilon )
    {
        displacementSumEq3.R = 0;
        displacementSumEq3.G = 0;
        displacementSumEq3.B = 0;
        displacementSumEq3.A = 0;

        //----------
        R_eq4.R = 0;
        R_eq4.G = 0;
        R_eq4.B = 0;
        R_eq4.A = 0;

        R_eq3.R = 0;
        R_eq3.G = 0;
        R_eq3.B = 0;
        R_eq3.A = 0;
        //----------

        displacementSumEq4.R = 0;
        displacementSumEq4.G = 0;
        displacementSumEq4.B = 0;
        displacementSumEq4.A = 0;
        return R_eq4;
    }

    //Section 3.5.1
    //During rendering in the fragment shader, for a given pixel (having texture coordinates x):
    //• For each channel aj : finalize the channel value computation by dividing the accumulated values by
    //sum(wi)
    //(resp. by sqrt(sum(wi2)) if using Eq. 4).
    //• Compute the texture value for the current pixel,
    //F (R(x)).
    //• Use the texture value as we would with a standard texture.

    //---------------- YU 2011 Equation 3 -----------------------
    //blending function, division part
    R_eq3.R = R_eq3.R/sumW;
    R_eq3.G = R_eq3.G/sumW;
    R_eq3.B = R_eq3.B/sumW;

    //---------------- YU 2011 Equation 4 -----------------------
    //blending function, division part
    float sqw = sqrtf(sumW2);

    if (outputCSV)
    {

        outfile.open("blending.csv", std::ios_base::app);
        vector<Pixel>::iterator itPixel;
        int i = 0;
        for(itPixel=colorsList.begin(); itPixel != colorsList.end(); itPixel++)
        {
            cout << "color "<<(*itPixel).R<<" "<<(*itPixel).G<<" "<<(*itPixel).B<<endl;
            //outfile << *itPixel<<",";
            cout << "w_i "<<w_i_list[i]<<endl;
            i++;
        }
        outfile <<sumW<<","<<sumW2<<","<< sqw << ","<<color_wi_sum.R << "," <<color_wi_sum.G << "," <<color_wi_sum.B << ","<<RM.R<<","<<RM.G<<","<<RM.B<< ",";
        cout << "sumW "<<sumW <<endl
                <<"sumW2 "<<sumW2 <<endl
                <<"sqw "<< sqw <<endl
                << "color_wi_sum "<<color_wi_sum.R << " " <<color_wi_sum.G << " " <<color_wi_sum.B <<endl
                << "RM "<<RM.R<<" "<<RM.G<<" "<<RM.B <<endl;
    }
    R_eq4.R = (color_wi_sum.R)/sqw    + RM.R;
    R_eq4.G = (color_wi_sum.G)/sqw    + RM.G;
    R_eq4.B = (color_wi_sum.B)/sqw    + RM.B;

    if (outputCSV)
    {
        outfile <<R_eq4.R << "," <<R_eq4.G << "," <<R_eq4.B<<endl;
        cout << R_eq4.R << " " <<R_eq4.G << " " <<R_eq4.B<<endl;
    }
    if (computeDisplacement)
    {
        displacementSumEq4.R = (displacementSumEq4.R)/sqw    + displaceMean.R;
        displacementSumEq4.G = (displacementSumEq4.G)/sqw    + displaceMean.G;
        displacementSumEq4.B = (displacementSumEq4.B)/sqw    + displaceMean.B;

        Clamp(displacementSumEq4);
    }

    float Af = sqw;
    if (Af > 1)
        Af = 1;
    //If we don't use the alpha, we are having black spot where sumW is 0
    //R_eq4.A = Af;
    //R1.A = Af;
    //we need to invertigate why we need to check that:
    Clamp(R_eq4);

    return R_eq4;
}




