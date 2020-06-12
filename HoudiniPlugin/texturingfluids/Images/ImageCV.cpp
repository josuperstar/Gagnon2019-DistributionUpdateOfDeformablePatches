#include "ImageCV.h"
#include <opencv2/imgproc/imgproc.hpp>

ImageCV::~ImageCV()
{
    //cout << "Destroying Image CV "<<this->m_id<< " "<< this->filename<<endl;
    if (this->image.dims != 0)
        this->image.release();
}

    bool ImageCV::OpenImage(const string& file_name, long id)
    {
        //cout << "Open image "<<id<< " "<< file_name <<" : ";
        cv::Mat temp = cv::imread(file_name, cv::IMREAD_UNCHANGED);
        if (temp.empty())
        {
            cout << " file "<<file_name<<" does not exist"<<endl;
            return false;
        }

        //cv::imshow("test",temp);
        this->image = temp;
        //cv::waitKey();
        cv::cvtColor(temp,this->image,cv::COLOR_BGR2BGRA);
        //cv::cvtColor(temp,this->image, cv::IMREAD_UNCHANGED);

        //m_width = image.size()[1];
        this->m_width = this->image.size().width;
        this->m_height = this->image.size().height;
        this->m_depth = 1;

        this->filename = file_name;

        this->m_id = id;
        this->isValid = true;
        //cout << " Ok"<<endl;
        return true;
    }

    bool ImageCV::OpenPNM(const string& file_name, long id)
    {
        //cout << "Open image "<< file_name <<" : ";
        cv::Mat temp = cv::imread(file_name);
        if (temp.empty())
        {
            cout << " file "<<file_name<<" does not exist"<<endl;
            return false;
        }

        //cv::imshow("test",temp);
        this->image = temp;
        //cv::waitKey();
        cv::cvtColor(temp,this->image,cv::COLOR_GRAY2RGBA);
        //m_width = image.size()[1];
        this->m_width = this->image.size().width;
        this->m_height = this->image.size().height;

        this->filename = file_name;
        this->m_depth = 1;
        this->m_id = id;
        this->isValid = true;
       // cout << " Ok"<<endl;
        return true;
    }

    void ImageCV::CreateImage(int dimension,vector<Pixel> pixels, bool useCenterPixel)
    {
        //int size = sqrt(pixels.size());

        int centerPosition = (dimension-1)/2;

        CreateImage(dimension, dimension, -1);
        for(int j=0;j<dimension;j++)
        {
            for(int i=0;i<dimension;i++)
            {
                if (!useCenterPixel && i == centerPosition && j == centerPosition && dimension != 1)
                    continue;

                this->SetColor(i,j,0,pixels[i+dimension*j]);
            }
        }


    }

    void ImageCV::OpenImageOrCreateImage(const string& file_name, int Width, int Height, long id)
    {
        //cout << "OpenImageOrCreateImage "<<id<<endl;
        //cv::Mat temp = cv::imread(file_name,CV_LOAD_IMAGE_UNCHANGED);
        cv::Mat temp = cv::imread(file_name);
        if (temp.empty())
        {
            CreateImage(Width, Height, id);
        }
        else
        {
            //cout << "open image "<<file_name<<endl;
            this->image = temp;
            cv::cvtColor(temp,this->image,cv::COLOR_BGR2RGB);
            m_width = this->image.size().width;
            m_height = this->image.size().height;
            this->m_depth = 1;
            //cout << "create image with "<<m_width<< " by "<<m_height<<endl;
            this->filename = file_name;
        }
        this->filename = file_name;
        this->isValid = true;
    }
		
    void ImageCV::CreateImage(int Width, int Height, long id)
    {
        //cout << "Create blank image "<<id<<" "<<Width<<"x"<<Height<<endl;
        //image.create(cv::Size(Width, Height),	CV_16UC4);
        //image.create(cv::Size(Width, Height), CV_8UC3);
        //cv::Mat temp;
        this->image.release();
        this->image.create(cv::Size(Width, Height), CV_8UC4);
        this->image = cv::Scalar(0,0,0,0);
        //cv::cvtColor(temp,image,CV_BGR2BGRA);
        this->m_width = Width;
        this->m_height = Height;
        this->m_depth = 1;
        this->m_id = id;
        this->isValid = true;

    }
    void ImageCV::CreateImage(int Width, int Height, const string& file_name, long id)
    {

        CreateImage(Width,Height, id);
        this->filename = file_name;

    }
    void ImageCV::SaveImage()
    {
        //cv::imshow("test", image);
        //cv::waitKey();

        if (this->image.rows == 0 || this->image.cols == 0)
        {
            cout << "can't save image, invalid size "<<endl;
            return;
        }

        vector<int> compression_params;
        //compression_params.push_back(CV_IMWRITE_PNG_COMPRESSION);
        //compression_params.push_back(9);
        //cv::imwrite(this->filename,this->image,compression_params);
        cv::imwrite(this->filename,this->image);
        //this->image.release();
        //growRegions(image,image,2);

        //cv::imwrite(file_name+"dilated.tiff",image);

    }

    void ImageCV::SaveImageAs(const string& file_name)
    {
        if (this->image.rows == 0 || this->image.cols == 0)
        {
            cout << "error in image "<<endl;
            return;
        }
        cv::imwrite(file_name,this->image);
    }

    void ImageCV::Clear()
    {
        //cout << "clearing image "<<this->m_id<<" "<<this->filename<<endl;
        this->image.release();
    }

    void ImageCV::GetColor(int x, int y,int z, Pixel& color)
    {

        if (this->image.empty())
        {
            cout << "no image "<<this->filename<<endl;
            return;
        }

        AdjustCoordinates(x,y);

        if (x < 0 || y < 0 || x >= this->m_width || y >=this->m_height || y >= this->image.size().height|| x >= this->image.size().width)
        {
            cout << "This is wrong !!!"<<endl;
            return;
        }


        float r = (float)this->image.at<cv::Vec4b>(y,x)[2];
        float g = (float)this->image.at<cv::Vec4b>(y,x)[1];
        float b = (float)this->image.at<cv::Vec4b>(y,x)[0];
        float a = (float)this->image.at<cv::Vec4b>(y,x)[3];
        //float a = 1;
        r /= 255.0f;
        g /= 255.0f;
        b /= 255.0f;
        a /= 255.0f;

        color.R = r;
        color.G = g;
        color.B = b;
        color.A = a;

        color.pixelPosition = Vec3f(x,y,z);

        //cout << "Get color " <<x<<"x"<<y<<" "<<color->R<<" "<<color->G<<" "<<color->B<<" "<<color->A<<" "<<endl;
    }

    void ImageCV::SetColor(int x, int y, int z, Pixel &color)
    {

        if (this == 0x0)
            return;
        if (this->image.cols <= 0 || this->image.rows <=0)
        {
            //cout << "problem with the image"<<endl;
            return;
        }

        AdjustCoordinates(x,y);

        if (x < 0 || y < 0 || x >= this->m_width || y >=this->m_height || y >= this->image.size().height|| x >= this->image.size().width)
        {
            cout << "This is wrong !!!"<<endl;
            return;
        }

        //if (x == 1642 && y == 1485)

        float r = color.R * 255.0f;
        float g = color.G * 255.0f;
        float b = color.B * 255.0f;
        float a = color.A * 255.0f;

        //cout << "Set color " <<x<<"x"<<y<<" "<<color->R<<" "<<color->G<<" "<<color->B<<" "<<color->A<<" "<<endl;
        /*
        image.at<cv::Vec4b>(x,y)[2] = color->R;
        image.at<cv::Vec4b>(x,y)[1] = color->G;
        image.at<cv::Vec4b>(x,y)[0] = color->B;
        image.at<cv::Vec4b>(x,y)[3] = color->A;
        */

        //cout << "color "<< color->R <<" "<<color->G<<" "<<color->B<<" "<<color->A<<endl;

        this->image.at<cv::Vec4b>(y,x)[2] = (uchar)r;
        this->image.at<cv::Vec4b>(y,x)[1] = (uchar)g;
        this->image.at<cv::Vec4b>(y,x)[0] = (uchar)b;
        this->image.at<cv::Vec4b>(y,x)[3] = (uchar)a;
        //image.at<cv::Vec4b>(y,x)[3] = 255;

    }

    void ImageCV::AdjustCoordinates(int &x, int &y)
    {
        while (x < 0)
            x = this->image.size().width + x;
        //if (x < 0)
        //    x = 0;

        while (y < 0)
            y = this->image.size().height + y;
        //if (y < 0)
        //    y = 0;

        while (x >= this->image.size().width )
            x -= this->image.size().width;
        //if (x >= this->m_width)
        //    x = this->m_width-1;

        while (y >= this->image.size().height )
            y -= this->image.size().height;
        //if (y >= this->m_height )
        //    y = this->m_height-1;
    }

    void ImageCV::growRegions(cv::Mat& input, cv::Mat& output, int n)
	{

        //===========================================================
       // Isolate each image channel
       cv::Mat cn[input.channels()];
       cv::split(input, cn);

       // Locate the pixels to add
       cv::Mat dilatedAlpha;
       cv::Mat strel = cv::Mat::ones(cv::Size(2*n+1,2*n+1),cn[3].type()); // Define a 3x3 structuring element
       cv::dilate(cn[3], dilatedAlpha, strel);
       cv::Mat pixelsToAdd = (dilatedAlpha-cn[3]) > 127;

       // Generate the pixels to add
       cv::Mat dilatedPixels;
       //cv::morphologyEx(input, dilatedPixels, CV_MOP_DILATE, strel);
        cv::morphologyEx(input, dilatedPixels, 1, strel);

       // Combine the results
       input.copyTo(output);
       dilatedPixels.copyTo(output, pixelsToAdd);
       //==============================================================

	}

    ImageCV ImageCV::DownSample(long id)
    {

        //long id = -1;

        cv::Mat dst;
        ImageCV imageDst;
        //=========GAUSSIAN PYRAMID LEVEL 2 ============
        pyrDown(this->image, dst, cv::Size( this->image.cols/2, this->image.rows/2 ));

        imageDst.CreateImage(this->image.cols/2,this->image.rows/2, id);
        imageDst.image = dst;
        return imageDst;

    }

    ImageCV ImageCV::UpSample()
    {

        long id = -1;

        cv::Mat dst;
        ImageCV imageDst;
        pyrUp(this->image, dst, cv::Size( this->image.cols*2, this->image.rows*2 ));

        imageDst.CreateImage(this->image.cols*2,this->image.rows*2, id);
        imageDst.image = dst;
        return imageDst;
    }

    vector<ImageCV> ImageCV::GaussianPyramid(int numberOfLevel)
    {

        vector<ImageCV> results;
        int level = numberOfLevel;
        vector<cv::Mat> gaussianPyramid;
        //vector<cv::Mat> laplacianPyramid;
        cv::Mat temp1, temp2, temp3;
        //cv::Mat Lap;
        this->image.copyTo(temp1);

        //level 0;
        ImageCV imageDist;
        imageDist.CreateImage(temp1.cols,temp1.rows,-1);
        imageDist.image = temp1;
        results.push_back(imageDist);

        for(int i=0; i<(level-1);i++)
        {
            pyrDown(temp1,temp2);
            //pyrUp(temp2,temp3,temp1.size());
            //Lap = temp1-temp3;
            gaussianPyramid.push_back(temp2);
            //laplacianPyramid.push_back(Lap);
            temp1 = temp2;

            //-----------------------------------

            imageDist.CreateImage(temp2.cols,temp2.rows,-1);
            imageDist.image = temp2;
            results.push_back(imageDist);

        }
    return results;
    }

    void ImageCV::Blur()
    {

        cv::Mat dst;
        //cout << "Blur"<<endl;
        long currentId = this->m_id;
        ImageCV temp = this->DownSample(currentId);
        cv::Mat src = temp.image;
        //cv::Mat src = this->image;
        cv::resize(src,dst,cv::Size(),2,2,cv::INTER_LINEAR);
        //cv::resize(src,dst,Size(),2,2,cv::INTER_NEAREST);

        /*
        cv::Mat dst;
        int MAX_KERNEL_LENGTH = 3;
        //int DELAY_BLUR = 100;
        for ( int i = 1; i < MAX_KERNEL_LENGTH; i = i + 2 )
        {
            GaussianBlur( this->image, dst, Size( i, i ), 0, 0 );
            //if( display_dst( DELAY_BLUR ) != 0 ) { return 0; }
        }
        */
        this->image = dst;
    }

    Pixel ImageCV::GetRandomColor(long seed)
    {
        Vec3d pixelPosition;

        srand(seed);
        int randomx = rand() % m_width;
        srand(seed+10);
        int randomy = rand() % m_height;
        pixelPosition.x = randomx;
        pixelPosition.y = randomy;
        Pixel c;

        this->GetColor(pixelPosition.x,pixelPosition.y,0,c);

        return c;
    }

    Pixel ImageCV::MeanValue()
    {
        cv::Scalar s = cv::mean(this->image);
        Pixel r;
        r.R = (s[2])/255;
        r.G = (s[1])/255;
        r.B = (s[0])/255;
        return r;

    }
