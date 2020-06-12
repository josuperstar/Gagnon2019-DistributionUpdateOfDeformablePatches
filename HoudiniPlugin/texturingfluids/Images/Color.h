
#ifndef __PIXEL_H__
#define __PIXEL_H__

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <Math/Vec3.h>
#include <cmath>
using namespace std;
using namespace TexturingFluids;

class Pixel
{
	public:
        float R, G, B, A;
        Vec3f pixelPosition;
        Vec3f meshPosition;


        bool valide;


        Pixel(float r = 1.0f, float g = 1.0f, float b = 1.0f, float a = 1.0f)
		{
            this->R = r;
            this->G = g;
            this->B = b;
            this->A = a;
            this->pixelPosition = Vec3f(0,0,0);
            this->meshPosition = Vec3f(0,0,0);
            this->valide = true;
		}

		uint32_t
		ToUInt32() const
		{
            uint32_t r = (uint32_t)(this->R * 255.0f);
            uint32_t g = (uint32_t)(this->G * 255.0f);
            uint32_t b = (uint32_t)(this->B * 255.0f);
            uint32_t a = (uint32_t)(this->A * 255.0f);

			return (a << 24) | (r << 16) | (g << 8) | b;
		}

        /*Pixel
        operator + (const Pixel &c) const
		{
            return Pixel(this->R + c.R, this->G + c.G, this->B + c.B, this->A + c.A);
        }*/
        Pixel
        operator + (const Pixel &c) const
        {
            return Pixel(R + c.R, G + c.G, B + c.B, A + c.A);
        }


        Pixel
        operator - (const Pixel &c) const
		{
            return Pixel(this->R - c.R, this->G - c.G, this->B - c.B, this->A - c.A);
		}

        Pixel
		operator * (float f) const
		{
            return Pixel(this->R * f, this->G * f, this->B * f, this->A * f);
		}

        Pixel
        operator * (Pixel p) const
        {
            return Pixel(this->R * p.R, this->G * p.G, this->B * p.B, this->A * p.A);
        }

        float Diff(const Pixel &c)
        {
            return fabs(this->R-c.R)+fabs(this->G-c.G)+fabs(this->B-c.B);
        }

		void Print()
		{
            std::cout << "["<<this->R<<" , "<<this->G<<" , "<<this->B<<"] mesh position "<<this->meshPosition<<" pixel position "<<this->pixelPosition<< std::endl;
		}


};


#endif /* __COLOR_H__ */
