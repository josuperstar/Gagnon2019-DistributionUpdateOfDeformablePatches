#ifndef __TBB_Atlas_h__
#define __TBB_Atlas_h__

#include "AtlasYu2011.h"

using namespace TexturingFluids;

struct Yu2011_executor
{
  Yu2011_executor(HoudiniAtlas &rasterizer, int w, int h,
           ParametersDeformablePatches params) : _rasterizer(rasterizer),
            _w(w), _h(h), _params(params)
  {
  }

  void operator()(const tbb::blocked_range<size_t>& r) const
  {
    //cout << "TBB Atlas"<<endl;
    for (size_t i=r.begin();i!=r.end();++i)
    {
      _rasterizer.RasterizePrimitive(GA_Offset(i), _w, _h, _params);
    }
  }

  HoudiniAtlas& _rasterizer;
  int _w;
  int _h;
  ParametersDeformablePatches _params;

};

#endif
