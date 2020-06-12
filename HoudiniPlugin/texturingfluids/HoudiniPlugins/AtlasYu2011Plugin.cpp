
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
#include <Core/HoudiniUtils.h>
#include <Strategies/StrategyPatchSurfaceSynthesis.h>
#include "AtlasYu2011Plugin.h"

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
    PRM_Name("TextureAtlasWidth",	"Texture Atlas Width"),
    PRM_Name("TextureAtlasHeight",	"Texture Atlas Height"),
    PRM_Name("TextureExemplar1",	"Texture Exemplar 1"),
    PRM_Name("TextureExemplarMask1",	"Texture Exemplar Mask 1"),
    PRM_Name("DisplacementMap1",	"Displacement Map 1"),
    PRM_Name("ComputeAtlas",	"Compute Atlas"),                   //5
    PRM_Name("TrackersFilename",	"Trackers Filename"),
    PRM_Name("RenderColoredPatches","Render Colored Patches"),
    PRM_Name("UseDeformableGrids","Use Deformable Grids"),
    PRM_Name("OutputName","Output Name"),
    PRM_Name("PoissonDiskRadius",	"Poisson Disk Radius"),
    PRM_Name("UVScaling",	"UV Scaling"),
    PRM_Name("PatchScaling",	"Patch Scaling"),
};

static PRM_Default PatchScalingDefault(1.0f);

PRM_Template
AtlasYu2011Plugin::myTemplateList[] =
{
    PRM_Template(PRM_TOGGLE, 1, &names[5]),
    PRM_Template(PRM_INT, 1, &names[0]),
    PRM_Template(PRM_INT, 1, &names[1]),
    PRM_Template(PRM_PICFILE_E, 1, &names[2]),
    PRM_Template(PRM_PICFILE_E, 1, &names[3]),
    PRM_Template(PRM_PICFILE_E, 1, &names[4]),
    PRM_Template(PRM_GEOFILE, 1, &names[6]),
    PRM_Template(PRM_TOGGLE, 1, &names[7]),
    PRM_Template(PRM_TOGGLE, 1, &names[8]),
    PRM_Template(PRM_STRING, 1, &names[9]),
    PRM_Template(PRM_FLT, 1, &names[10]),
    PRM_Template(PRM_FLT, 1, &names[11]),
    PRM_Template(PRM_FLT, 1, &names[12], &PatchScalingDefault),
    PRM_Template(),
};


OP_Node *
AtlasYu2011Plugin::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new AtlasYu2011Plugin(net, name, op);
}

AtlasYu2011Plugin::AtlasYu2011Plugin(OP_Network *net, const char *name, OP_Operator *op)
	: SOP_Node(net, name, op), myGroup(0)
{
    // Make sure to flag that we can supply a guide geometry
    mySopFlags.setNeedGuide1(1);
}

AtlasYu2011Plugin::~AtlasYu2011Plugin()
{
    cout << "Destroying DeformablePatches"<<endl;
    //this->interface.~UnitTestInterface();
}

OP_ERROR
AtlasYu2011Plugin::cookInputGroups(OP_Context &context, int alone)
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
AtlasYu2011Plugin::cookMySop(OP_Context &context)
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
    //int numberOfGaussianLevel = 2;

    cout << "======================== Gagnon2019DistributionUpdateOfDeformablePatches Atlas ============================="<<endl;

    string baseVariable = "REZ_GAGNON2019_DISTRIBUTION_UPDATE_DEFORMABLE_PATCHES_BASE";
    char* pPath;
    pPath = getenv (baseVariable.c_str());
    if (pPath!=NULL)
    cout << "version "<<pPath<<endl;
    //int startFrame = StartFrame();
    //int startNumber = 0;
    ParametersDeformablePatches params;

    params.frame = frame;
    params.poissondiskradius = PoissonDiskRadius();
    params.atlasHeight = TextureAtlasHeight();
    params.atlasWidth = TextureAtlasWidth();

    params.useDeformableGrids = UseDeformableGrids();
    params.coloredPatches = RenderColoredPatches();
    params.UVScaling = UVScaling(now);
    params.PatchScaling = PatchScaling(now);

    params.NumberOfTextureSampleFrame = 1;

    if (params.atlasHeight <= 0)
    {
        params.atlasHeight = 100;
    }
    if (params.atlasWidth <= 0)
    {
        params.atlasWidth = 100;
    }

    cout << "frame :"<<frame<<endl;
    cout << "disk radius :"<<params.poissondiskradius<<endl;
    cout << "atlasHeight :"<<params.atlasHeight<<endl;
    cout << "atlasWidth :"<<params.atlasWidth<<endl;
    cout << "UVScaling :"<<params.UVScaling<<endl;
    cout << "PatchScaling :"<<params.PatchScaling<<endl;

    TrackersFilename(trackersFilename,now);
    params.trackersFilename = trackersFilename;
    /*
    DeformableGridsFilename(deformableGridsFilename,now);
    params.deformableGridsFilename = deformableGridsFilename;
    */
    TextureExemplar1(textureExemplar1Name,now);
    params.textureExemplar1Name = textureExemplar1Name;

    TextureExemplarMask1(textureExemplar1MaskName,now);
    params.textureExemplar1MaskName = textureExemplar1MaskName;

    DisplacementMap1(displacementMap1Name,now);
    params.displacementMap1Name = displacementMap1Name;

    OutputName(outputName,now);
    params.outputName = outputName;

    const GU_Detail * surface = inputGeo(1);
    surface = inputGeo(1);

    GU_Detail *surfaceCopy = new GU_Detail();
    surfaceCopy->clearAndDestroy();
    surfaceCopy->copy(*surface);

    const GU_Detail *trackersGdp = inputGeo(2);
    GU_Detail *trackersCopy = new GU_Detail();
    trackersCopy->clearAndDestroy();
    trackersCopy->copy(*trackersGdp);

    AtlasInterface interface;
    bool synthesised = interface.Synthesis(gdp, surfaceCopy,trackersCopy, params);
    if (synthesised)
        cout << "was able to synthesis the atlas"<<endl;
    else
        cout << "was not able to synthesis the atlas"<<endl;


    unlockInputs();
    resetLocalVarRefs();

    cout << "==================================================================="<<endl;
    return error();
}



const char *
AtlasYu2011Plugin::inputLabel(unsigned) const
{
    return "Surface Deformable Patches";
}
