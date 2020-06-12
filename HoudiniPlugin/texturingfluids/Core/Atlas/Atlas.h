#ifndef __Atlas_h__
#define __Atlas_h__

#include <string>
#include <vector>

namespace TexturingFluids {


class Atlas
{

public:

    virtual bool BuildAtlas(int w, int h, int life) = 0;


    void SetFilename(std::string data){outputFilename = data;}
    std::string GetFilename() {return outputFilename;}

protected :
    std::string outputFilename;


	
};
} // End HDK_Sample namespace

#endif
