#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include "Vec3.h"
#include "math.h"

using namespace std;


namespace TexturingFluids {

class Triangle{
	public :

    bool static SameSide(TexturingFluids::Vec3d p1,TexturingFluids::Vec3d p2,TexturingFluids::Vec3d a, TexturingFluids::Vec3d b);
    bool static PointInTriangle2d(TexturingFluids::Vec3d p, TexturingFluids::Vec3d a,TexturingFluids::Vec3d b,TexturingFluids::Vec3d c);
    TexturingFluids::Vec3d static CreateVector( TexturingFluids::Vec3d a,TexturingFluids::Vec3d b);
    static double TriangleAir(TexturingFluids::Vec3d a,TexturingFluids::Vec3d b,TexturingFluids::Vec3d c);
    bool static IsPointInTriangle(TexturingFluids::Vec3d  p, TexturingFluids::Vec3d a,TexturingFluids::Vec3d b,TexturingFluids::Vec3d c);
    bool static IsPointInTriangle2(TexturingFluids::Vec3d  p, TexturingFluids::Vec3d a,TexturingFluids::Vec3d b,TexturingFluids::Vec3d c);
    bool static IsPointInTriangle2D(const TexturingFluids::Vec3d& p, const TexturingFluids::Vec3d& a, const TexturingFluids::Vec3d& b, const TexturingFluids::Vec3d& c);
    bool static IsPointInTriangleRange(TexturingFluids::Vec3d  p, TexturingFluids::Vec3d a,TexturingFluids::Vec3d b,TexturingFluids::Vec3d c, float s);
    static TexturingFluids::Vec3d ProjectPointOnTriangle(TexturingFluids::Vec3d p,  TexturingFluids::Vec3d a,TexturingFluids::Vec3d b,TexturingFluids::Vec3d c);
    TexturingFluids::Vec3d static GetUVBarycenter(TexturingFluids::Vec3d a,TexturingFluids::Vec3d b, TexturingFluids::Vec3d c,TexturingFluids::Vec3d A,TexturingFluids::Vec3d B, TexturingFluids::Vec3d C, TexturingFluids::Vec3d P);
    static TexturingFluids::Vec3d GetBarycentricPosition(TexturingFluids::Vec3d A,TexturingFluids::Vec3d B, TexturingFluids::Vec3d C, TexturingFluids::Vec3d a, TexturingFluids::Vec3d b, TexturingFluids::Vec3d c, TexturingFluids::Vec3d position);
	inline
        static double squaredDistToEdge(const TexturingFluids::Vec3d& p,
                                 const TexturingFluids::Vec3d& v0,
                                 const TexturingFluids::Vec3d& edgeDirection,
								 double			  edgeLength)
		{
			// Project point onto the edge
            TexturingFluids::Vec3d v0p(p);
			v0p -= v0;

			double d = v0p.x*edgeDirection.x + v0p.y*edgeDirection.y + v0p.z*edgeDirection.z;

			if (d<0.0) d = 0.0;
			if (d>edgeLength) d = edgeLength;

            TexturingFluids::Vec3d projP(edgeDirection);
			projP *= d;

			// Compute the squared distance
            TexturingFluids::Vec3d delta(projP);
			delta -= v0p;

			return delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
		}
    void static BoundingBox(TexturingFluids::Vec3d a,TexturingFluids::Vec3d b,TexturingFluids::Vec3d c, TexturingFluids::Vec3d &min,TexturingFluids::Vec3d &max);
    void static BoundingBox2D(TexturingFluids::Vec3d a,TexturingFluids::Vec3d b,TexturingFluids::Vec3d c, TexturingFluids::Vec3d &min,TexturingFluids::Vec3d &max);
};

}

#endif
