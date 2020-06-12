
#include <vector>
#include <algorithm>
#include <SYS/SYS_Math.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Matrix4.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_SpareData.h>
#include <SOP/SOP_Guide.h>
//#include "DeformablePatches.h"
#include <UT/UT_DSOVersion.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <GA/GA_ElementWrangler.h>
#include <algorithm>
#include <ctime>
#include <Core/HoudiniUtils.h>
#include <Strategies/StrategyPatchSurfaceSynthesis.h>
#include "Yu2011Plugin.h"


#include <stdlib.h> /* getenv */


#include <omp.h>

//#include "vector.h"
#include <Math/Vec3.h>

#include <iostream>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/convex_hull_2.h>

#define DEBUG 0

using namespace std;

using namespace TexturingFluids;


//===================================================================================================================
//===================================================================================================================
//===================================================================================================================



static PRM_Name        names[] = {

    PRM_Name("StartFrame",	"Start Frame"),
    PRM_Name("UpdateDistribution",	"Update Distribution"),
    PRM_Name("ComputeDistortion",	"ComputeDistortion"),
    PRM_Name("UseDeformableGrids",	"Use Deformable Grids for Atlas"),
    PRM_Name("MinimumDistanceProjection",	"Minimum Distance Projection"),
    PRM_Name("PoissonDiskRadius",	"Poisson Disk Radius"),//5
    PRM_Name("FadingTau",	"Fading Tau"),
    PRM_Name("Yu2011DMax",	"Yu 2011 DMax"),
    PRM_Name("QvMin",	"Yu 2011 QvMin"),
    PRM_Name("TestPatch", "Test Patch"),
    PRM_Name("PatchNumber",	"PatchNumber"),//10
    PRM_Name("TrackersFilename",	"Trackers Filename"),
    PRM_Name("DeformableGridsFilename",	"Deformable Grids Filename"),
    PRM_Name("PatchAngleNormalThreshold",	"Patch Angle Normal Threshold"),
    PRM_Name("PoissonAngleNormalThreshold",	"Poisson Angle Normal Threshold"),
    PRM_Name("UVScaling",	"UV Scaling"),
    PRM_Name("Yu2011Beta",	"Yu2011 Beta"),
    PRM_Name("CellSize",	"Cell Size"),
    PRM_Name("UseDynamicFading",	"Use Dynamic Fading"),
    PRM_Name("FadingIn",	"Use Fading In"),
};

static PRM_Default StartFrameDefault(1);
static PRM_Default PoissonDiskRadiusDefault(1.0f);
static PRM_Default MinimumDistanceProjectionDefault(0.01f);
static PRM_Default FadingTauDefault(48);
static PRM_Default Yu2011DMaxDefault(2.0f);
static PRM_Default QvMinDefault(0.5f);
static PRM_Default AngleNormalThresholdDefault(0.5f);
static PRM_Default PoissonAngleNormalThresholdDefault(0.9f);
static PRM_Default UVScalingDefault(2.0f);
static PRM_Default Yu2011BetaDefault(0.6f);
static PRM_Default CellSizeDefault(0.1f);
static PRM_Default UseDynamicFadingDefault(1);
static PRM_Default FadingInDefault(1);

PRM_Template
LagrangianTextureAdvectionPlugin::myTemplateList[] = {
    PRM_Template(PRM_FLT, 1, &names[0], &StartFrameDefault),
    PRM_Template(PRM_TOGGLE, 1, &names[1]),
    PRM_Template(PRM_TOGGLE, 1, &names[2]),
    PRM_Template(PRM_TOGGLE, 1, &names[3]),
    PRM_Template(PRM_FLT, 1, &names[4], &MinimumDistanceProjectionDefault),
    PRM_Template(PRM_FLT, 1, &names[5],&PoissonDiskRadiusDefault),
    PRM_Template(PRM_FLT, 1, &names[6], &FadingTauDefault),
    PRM_Template(PRM_FLT, 1, &names[7]),
    PRM_Template(PRM_FLT, 1, &names[8]),
    PRM_Template(PRM_TOGGLE, 1, &names[9]),
    PRM_Template(PRM_INT, 1, &names[10]),
    PRM_Template(PRM_GEOFILE, 1, &names[11]),
    PRM_Template(PRM_GEOFILE, 1, &names[12]),
    PRM_Template(PRM_FLT, 1, &names[13], &AngleNormalThresholdDefault),
    PRM_Template(PRM_FLT, 1, &names[14], &PoissonAngleNormalThresholdDefault),
    PRM_Template(PRM_FLT, 1, &names[15], &UVScalingDefault),
    PRM_Template(PRM_FLT, 1, &names[16], &Yu2011BetaDefault),
    PRM_Template(PRM_FLT, 1, &names[17], &CellSizeDefault),
    PRM_Template(PRM_INT, 1, &names[19], &FadingInDefault),
    PRM_Template(PRM_TOGGLE, 1, &names[18]),
    PRM_Template(),

};


OP_Node *
LagrangianTextureAdvectionPlugin::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new LagrangianTextureAdvectionPlugin(net, name, op);
}

LagrangianTextureAdvectionPlugin::LagrangianTextureAdvectionPlugin(OP_Network *net, const char *name, OP_Operator *op)
    : SOP_Node(net, name, op), myGroup(0)
{
    // Make sure to flag that we can supply a guide geometry
    mySopFlags.setNeedGuide1(1);
}

LagrangianTextureAdvectionPlugin::~LagrangianTextureAdvectionPlugin()
{
    cout << "Destroying DeformablePatches"<<endl;
    //this->interface.~UnitTestInterface();
}

OP_ERROR
LagrangianTextureAdvectionPlugin::cookInputGroups(OP_Context &context, int alone)
{
    // If we are called by the handle, then "alone" equals 1.  In that
    // case, we have to lock the inputs oursevles, and unlock them
    // before exiting this method.
    if (alone) if (lockInputs(context) >= UT_ERROR_ABORT) return error();

    UT_String	 grp_name;

    // The "gdp" variable is only available if we are called from the SOP
    // itself.  So, if we are called by a handle, we have to get the
    // geometry oursevles.
    GU_Detail	*pgdp = alone ? (GU_Detail *)inputGeo(0, context) : gdp;

    myGroup = 0;
    primGroup = 0;

    getGroups(grp_name);		// Get the group string.

    // If the group string is not null, then we try to parse the group.
    if (grp_name.isstring())
    {
    myGroup = parseEdgeGroups((const char *)grp_name, pgdp);

    // If the group is not valid, then the group string is invalid
    // as well.  Thus, we add an error to this SOP.
    if (!myGroup)
    {
        addError(SOP_ERR_BADGROUP, grp_name);
    }
    else if (!alone)
    {
        // If the parsed group is valid, then we want to highlight
        // only the group.  The second argument of "1" means that
        // we want the selection to have the same type as our group.
        select(*const_cast<GA_EdgeGroup*>(myGroup), 1);
    }
    }
    else if (!alone)
    {
    // If no group string is specified, then we operate on the entire
    // geometry, so we highlight every point for this SOP.
    select(GU_SPoint);
    }

    // This is where we notify our handles (if any) if the inputs have changed.
    //checkInputChanged(0, -1, myDetailGroupPair, pgdp, myGroup);

    // If we are called by the handles, then we have to unlock our inputs.
    if (alone)
    {
    destroyAdhocGroups();
    unlockInputs();
    }

    return error();
}


OP_ERROR
LagrangianTextureAdvectionPlugin::cookMySop(OP_Context &context)
{
    // Before we do anything, we must lock our inputs.  Before returning,
    //	we have to make sure that the inputs get unlocked.
    if (lockInputs(context) >= UT_ERROR_ABORT)
    return error();

    duplicateSource(0, context);
    setVariableOrder(3, 2, 0, 1);
    setCurGdh(0, myGdpHandle);
    setupLocalVars();

    //float cellSize = 2;
    fpreal now = context.getTime();
    int frame = context.getFrame();



    string baseVariable = "REZ_GAGNON2019_DISTRIBUTION_UPDATE_DEFORMABLE_PATCHES_BASE";
    char* pPath;
    pPath = getenv (baseVariable.c_str());
    if (pPath!=NULL)
    cout << "version "<<pPath<<endl;

    int startFrame = StartFrame();
    int startNumber = 0;
    ParametersDeformablePatches params;
    params.frame = frame;
    params.startFrame = startFrame;
    params.startNumber = startNumber;
    params.poissondiskradius = PoissonDiskRadius();
    params.updateDistribution = UpdateDistribution();
    params.maximumProjectionDistance = MinimumDistanceProjection();
    params.computeDistortion = ComputeDistortion();
    params.fadingTau = FadingTau();
    params.Yu2011DMax = Yu2011DMax();
    params.QvMin = QvMin();
    params.testPatch = TestPatch();
    params.patchNumber = PatchNumber();
    params.useDeformableGrids = UseDeformableGrids();
    params.angleNormalThreshold = PatchAngleNormalThreshold();
    params.poissonAngleNormalThreshold = PoissonAngleNormalThreshold();
    params.UVScaling = UVScaling();
    params.Yu2011Beta = Yu2011Beta();
    params.CellSize = CellSize();
    TrackersFilename(trackersFilename,now);
    params.trackersFilename = trackersFilename;
    params.fadingIn = FadingIn();

    DeformableGridsFilename(deformableGridsFilename,now);
    params.deformableGridsFilename = deformableGridsFilename;

    params.useDynamicTau = UseDynamicFading();
    if (params.useDynamicTau)
    {
        cout << "======================== Gagnon2019DistributionUpdateOfDeformablePatches, frame  "<<frame<< "============================="<<endl;
    }
    else
    {
        cout << "======================== YU 2011 Lagrangian Texture, frame  "<<frame<< "============================="<<endl;
    }

    const GU_Detail * surface = inputGeo(1);
    surface = inputGeo(1);
    GU_Detail *surfaceCopy = new GU_Detail();
    surfaceCopy->clearAndDestroy();
    surfaceCopy->copy(*surface);

    const GU_Detail *trackersGdp = inputGeo(2);
    GU_Detail *trackersCopy = new GU_Detail();
    trackersCopy->clearAndDestroy();
    trackersCopy->copy(*trackersGdp);

    const GU_Detail *levelSetRef = inputGeo(3);
    GU_Detail *levelSet = new GU_Detail();
    levelSet->clearAndDestroy();
    levelSet->copy(*levelSetRef);

    const GU_Detail *surfaceLowResRef = inputGeo(4);
    GU_Detail *surfaceLowRes = new GU_Detail();
    surfaceLowRes->clearAndDestroy();
    surfaceLowRes->copy(*surfaceLowResRef);

    LagrangianTextureAdvection interface;
    //interface.Synthesis(gdp,const_cast<GU_Detail*>(surface), params);
    interface.Synthesis(gdp,surfaceCopy,trackersCopy,levelSet,surfaceLowRes, params);

    delete trackersCopy;
    delete surfaceCopy;
    delete levelSet;
    delete surfaceLowRes;

    unlockInputs();
    resetLocalVarRefs();
    cout << "=============================== END ===================================="<<endl;
    return error();
}

const char *
LagrangianTextureAdvectionPlugin::inputLabel(unsigned) const
{
    return "Surface Deformable Patches";
}
