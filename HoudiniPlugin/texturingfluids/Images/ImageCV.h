#ifndef __IMAGECV_h__
#define __IMAGECV_h__


#include "Color.h"
#include <opencv2/opencv.hpp>


#define USEVIRTUAL 1

//using namespace std;
//using namespace cv;
#if USEVIRTUAL
#include "Image.h"

class ImageCV : public Image
#else
class ImageCV
#endif
{
	public:
    ImageCV() {}
    virtual ~ImageCV();


    bool OpenImage(const std::string& file_name, long id);
    bool OpenPNM(const std::string& file_name, long id);
		
    void CreateImage(int Width, int Height, long id);
    void CreateImage(int Width, int Height,const std::string& file_name, long id);
    void CreateImage(int dimension,vector<Pixel> pixels, bool useCenterPixel);
    void OpenImageOrCreateImage(const std::string& file_name, int Width, int Height, long id);
    void SaveImage();
    void SaveImageAs(const std::string& file_name);
    void GetColor(int x, int y,int z, Pixel& color);
    void SetColor(int x, int y,int z, Pixel &color);
    static void growRegions(cv::Mat& input, cv::Mat& output, int n);

    ImageCV DownSample(long id);
    ImageCV UpSample();
    vector<ImageCV> GaussianPyramid(int numberOfLevel);
    Pixel GetRandomColor(long seed);
    Pixel MeanValue();

    void Blur();
    void Clear();
    cv::Mat image;
#if USEVIRTUAL
#else
    int GetWidth() {return m_width;}
    int GetHeight() {return m_height;}
    std::string GetFilename(){return filename;}
#endif
private :
    void AdjustCoordinates(int &x, int &y);

#if USEVIRTUAL

#else
protected :

    int m_width;
    int m_height;
    std::string filename;
#endif


		
};

#endif // __IMAGECV_h__
