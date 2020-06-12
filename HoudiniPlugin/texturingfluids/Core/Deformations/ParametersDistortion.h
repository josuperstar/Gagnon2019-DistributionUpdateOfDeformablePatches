#ifndef __ParametersDistortion__
#define __ParametersDistortion__
#include <vector>
#include <SOP/SOP_Node.h>
#include <GEO/GEO_PrimPart.h>
#include <GEO/GEO_PointTree.h>
//#include <GA_ElementGroup.h>
#include <vector>



namespace TexturingFluids {


struct ParametersDistortion
{
    float dilatationMin;
    float dilatationMax;
    float squeezeMin;
    float squeezeMax;
    float deletionLife;
    float distortionRatioThreshold;
    std::string alphaName;
    std::string temporalRemoveName;
    std::string randomThresholdDistortion;
    std::string initialVertexAngle;
    std::string distortionWeightName;
    std::string primLifeName;
    bool computeSqueeze = false;
    bool computeDilatation = false;
    bool computeShearing = false;
    float Yu2011DMax;
    float QvMin;
    int flagDistortedParticles = true;
};

}
#endif
