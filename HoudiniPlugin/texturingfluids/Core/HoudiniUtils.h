#ifndef __HOUDINIUTILS__
#define __HOUDINIUTILS__
#include <vector>
#include <SOP/SOP_Node.h>
#include <GEO/GEO_PrimPart.h>
#include <GEO/GEO_PointTree.h>
//#include <GA_ElementGroup.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>


namespace TexturingFluids {


class HoudiniUtils
{

public:

    static std::set<GA_Offset> GetNeighbors(GU_Detail *gdp, GA_Offset ppt);
    static std::set<GA_Offset> GetPrimitivesNeighbors(GU_Detail *gdp, GA_Offset ppt);

    static void AttributeTransfert(GU_Detail *gdp, GU_Detail *referenceGdp, GEO_PointTreeGAOffset &tree, GA_PointGroup *grp, float transfertRadius, GA_RWHandleV3 refAtt, GA_RWHandleV3 destAtt);
    static void AttributeTransfert(GU_Detail *gdp, GU_Detail *referenceGdp, GEO_PointTreeGAOffset &tree, GA_PointGroup *grp, float transfertRadius, GA_RWHandleV3 refAtt, GA_RWHandleV3 destAtt, float excludingDot);
    static void AttributeTransfert(GU_Detail *gdp, GEO_PointTreeGAOffset &tree, GA_PointGroup *fromGroup, GA_PointGroup *destGroup, float transfertRadius,
                                   GA_RWHandleV3 refAtt,
                                   GA_Attribute        *uvsAtt,
                                   const GA_AIFNumericArray *uvsArray);

    static void AttributeTransfert(GU_Detail *gdp, GEO_PointTreeGAOffset &tree, GA_PointGroup *fromGroup, GA_PointGroup *destGroup, float transfertRadius,
                                   GA_RWHandleF refAtt,
                                   GA_Attribute        *uvsAtt,
                                   const GA_AIFNumericArray *uvsArray);

    static void AttributeTransfert(GU_Detail *gdp, GU_Detail *referenceGdp, GEO_PointTreeGAOffset &tree, float transfertRadius, GA_RWHandleI refAtt, GA_RWHandleI destAtt);

	
    static UT_Vector3 GetBarycentricPosition(UT_Vector3 A,UT_Vector3 B, UT_Vector3 C, UT_Vector3 a, UT_Vector3 b, UT_Vector3 c, UT_Vector3 position);

    static std::string format_account_number(int acct_no) {
      std::ostringstream out;
      out << std::internal << std::setfill('0') << std::setw(4) << acct_no;
      return out.str();
    }

};
}

#endif
