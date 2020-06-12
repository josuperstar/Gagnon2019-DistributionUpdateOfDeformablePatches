#ifndef __IMAGE_h__
#define __IMAGE_h__

#include "Color.h"

//using namespace std;

class Image
{
	public:
    Image() {}
    virtual ~Image()
    {
        //cout << "Destroying Virtual Image"<<endl;
    }

    virtual bool OpenImage(const string& file_name, long id) {return false;}// =0;

    virtual void CreateImage(int Width, int Height, long id) {}//=0;
    virtual void CreateImage(int Width, int Height, int Depth, long id) {}//=0;
    virtual void CreateImage(int Width, int Height,const string& file_name, long id) {}//=0;
    virtual void OpenImageOrCreateImage(const string& file_name, int Width, int Height, long id) {}//=0;

    virtual void SaveImage() {}//=0;
    virtual void SaveImageAs(const string& file_name) {}//=0;

    virtual void GetColor(int x, int y, int z, Pixel& color) const { cout << "getcolor mother"<<endl;}//=0;
    virtual void SetColor(int x, int y, int z, Pixel &color) {}//=0;

    virtual void Clear() {}//=0;

    int GetWidth() const {return m_width;}
    int GetHeight() const {return m_height;}
    int GetDepth() const {return m_depth;}
    string GetFilename() const {return filename;}
    long GetId() const {return m_id;}
    bool IsValid(){return isValid;}

//private :
    //void AdjustCoordinates(int &x, int &y);

protected :

    int m_width;
    int m_height;
    int m_depth;
    string filename;
    long m_id;
    bool isValid = false;



		
};

#endif // __IMAGE_h__
