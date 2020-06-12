#ifndef __ParametersDeformablePatches_h__
#define __ParametersDeformablePatches_h__

using namespace std;

namespace TexturingFluids {
struct ParametersDeformablePatches
{
    float poissondiskradius;
    float innerCircleRadius;
    int startFrame;
    int startNumber;
    int frame;
    float tangentTrackerLenght;
    int updateDistribution;
    float alphaThreshold;
    float maximumProjectionDistance;
    int connectivityTest;
    int testPatch;
    int patchNumber;
    float velocityTransfertRadius;
    float alphaTransfertRadius;
    float uvTransfertRadius;
    float normalTransfertRadius;
    int deleteConcealedPatches;
    int computeDistortion;
    float distortionRatioThreshold; //the amount of distorted vertex allow per grid
    float angleNormalThreshold;     //the angle (dot product) with the normal of the tracker below which the vertex is considered as distorted
    float poissonAngleNormalThreshold;
    string trackersFilename;
    string deformableGridsFilename;
    bool computeAtlas;
    string textureExemplar1Name;
    string textureExemplar1MaskName;
    string displacementMap1Name;
    string outputName;
    int atlasWidth;
    int atlasHeight;
    bool useDeformableGrids;
    bool coloredPatches;
    float fadingTau; // Yu2011 Fading Tau, equation 6
    float Yu2011DMax; // Yu2011 delta max, equation 2
    float QvMin; //Yu2011 Quality Vertex Minimum, 0.5f by default
    float UVScaling; //scaling used for the patch definition of uv.
    float PatchScaling; //scaling used for the patch definition size.
    bool useDynamicTau;
    float Yu2011Beta;
    float CellSize;
    int NumberOfTextureSampleFrame = 1;
    int useTangeantTracker = 0;
    int fadingIn = 1;
};
}

#endif
