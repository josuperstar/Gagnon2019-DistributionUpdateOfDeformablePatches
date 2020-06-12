
#include "triangle.h"

using namespace TexturingFluids;

bool SameSide(Vec3d p1,Vec3d p2,Vec3d a, Vec3d b)
{
	/*
	UT_Vector3 cp1 = cross(b-a, p1-a);
	UT_Vector3 cp2 = cross(b-a, p2-a);
	if (dot(cp1, cp2) >= 0)
	 return true;
	else return false;
	*/

    Vec3d ab = b - a;
    Vec3d ap1 = p1 -a;
    Vec3d ap2 = p2 - a;
    Vec3d cp1 = ab.cross(ap1);
        Vec3d cp2 = ab.cross(ap2);
        if (cp1.dot(cp2) >= 0)
		return true;
   	else return false;
}

float TriangleAir(Vec3d a,Vec3d b,Vec3d c)
{
	/*
	UT_Vector3 x =(c-b);
	UT_Vector3 y =(c-a);
	UT_Vector3 z = cross(x,y);
	
    return 0.5* z.length();
	*/

    Vec3d x =(c-b);
    Vec3d y =(c-a);
    Vec3d z = x.cross(y);
	
    return 0.5* z.length();
}


Vec3d sous(Vec3d A, Vec3d B)
{
    Vec3d r;
	r.x = B.x - A.x;
	r.y = B.y - A.y;
	r.z = B.z - A.z;
	return r;
}

Vec3d add(Vec3d A, Vec3d B)
{
    Vec3d r;
	r.x = B.x + A.x;
	r.y = B.y + A.y;
	r.z = B.z + A.z;
	return r;
}


bool Triangle::PointInTriangle2d(Vec3d p, Vec3d a,Vec3d b,Vec3d c)
{
    if (SameSide(p,a, b,c) && SameSide(p,b, a,c) && SameSide(p,c, a,b))
	return true;
    else return false;
}

Vec3d Triangle::CreateVector( Vec3d a,Vec3d b)
{
    Vec3d result;
	result.x = b.x - a.x;
	result.y = b.y - a.y;
	result.z = b.z - a.z;
	return result;
}

double Triangle::TriangleAir(Vec3d a,Vec3d b,Vec3d c)
{
    double xx = c.x-b.x;
    double xy = c.y-b.y;
    double xz = c.z-b.z;

    double yx = c.x-a.x;
    double yy = c.y-a.y;
    double yz = c.z-a.z;

    double zx = xy*yz - xz*yy;
    double zy = xz*yx - xx*yz;
    double zz = xx*yy - xy*yx;

    return 0.5*sqrt(zx*zx + zy*zy + zz*zz);

    /*
    Vec3d x =(c-b);
    Vec3d y =(c-a);
    Vec3d z = x.cross(y);
	
    return 0.5* z.length();*/
}

bool Triangle::IsPointInTriangle2D(const Vec3d& p, const Vec3d& a, const Vec3d& b, const Vec3d& c)
{
	// Determine on wich side of each lines the point lies
	double sideAB = (p.x-a.x)*(b.y-a.y) - (p.y-a.y)*(b.x-a.x);
	double sideBC = (p.x-b.x)*(c.y-b.y) - (p.y-b.y)*(c.x-b.x);
	double sideCA = (p.x-c.x)*(a.y-c.y) - (p.y-c.y)*(a.x-c.x);

	// Are we always on the good side of the line?
	return (sideAB>=0.0 && sideBC>=0.0 && sideCA>=0.0);
}

bool Triangle::IsPointInTriangle(Vec3d  p, Vec3d a,Vec3d b,Vec3d c)
{
	/*
        // Compute vectors
        UT_Vector3 v0 = c - a;
        UT_Vector3 v1 = b - a;
        UT_Vector3 v2 = p - a;

        // Compute dot products
        float dot00 = dot(v0, v0);
        float dot01 = dot(v0, v1);
        float dot02 = dot(v0, v2);
        float dot11 = dot(v1, v1);
        float dot12 = dot(v1, v2);

        // Compute barycentric coordinates
        float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
        float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

        // Check if point is in triangle
        return (u >= 0) && (v >= 0) && (u + v <= 1);
	*/
	
    float epsilon = 0.00001;
    Vec3d v0 = c - a;
    Vec3d v1 = b - a;
    Vec3d v2 = p - a;


    // Compute dot products
    float dot00 = v0.dot(v0);
    float dot01 = v0.dot(v1);
    float dot02 = v0.dot(v2);
    float dot11 = v1.dot(v1);
    float dot12 = v1.dot(v2);

    // Compute barycentric coordinates
    float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Check if point is in triangle
    return (u >= 0) && (v >= 0) && (u + v <= 1);
}

bool Triangle::IsPointInTriangle2(Vec3d  p, Vec3d a,Vec3d b,Vec3d c)
{
	/*
		// Compute vectors
		UT_Vector3 v0 = t->c - t->a;
		UT_Vector3 v1 = t->b - t->a;
		UT_Vector3 v2 = p - t->a;

		// Compute dot products
		float dot00 = dot(v0, v0);
		float dot01 = dot(v0, v1);
		float dot02 = dot(v0, v2);
		float dot11 = dot(v1, v1);
		float dot12 = dot(v1, v2);

		// Compute barycentric coordinates
		float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
		float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
		float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

		// Check if point is in triangle
		return (u >= 0) && (v >= 0) && (u + v <= 1);
	*/

	float epsilon = 0.0001;
    Vec3d v0 = c - a;
    Vec3d v1 = b - a;
    Vec3d v2 = p - a;


    Vec3d e1 = b-a;
    Vec3d e2 = b-c;
    Vec3d e3 = a-c;

    Vec3d n1 = e1;
    Vec3d n2 = e2;
    Vec3d n3 = e3;

    n1.normalize();
    n2.normalize();
    n3.normalize();

    double d1 = squaredDistToEdge(p, a, n1, e1.length() );
    double d2 = squaredDistToEdge(p, b, n2, e2.length() );
    double d3 = squaredDistToEdge(p, c, n3, e3.length() );


	// Compute dot products
    float dot00 = v0.dot(v0);
    float dot01 = v0.dot(v1);
    float dot02 = v0.dot(v2);
    float dot11 = v1.dot(v1);
    float dot12 = v1.dot(v2);

	// Compute barycentric coordinates
	float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point is in triangle
	return (((u >= 0) && (v >= 0) && (u + v <= 1)) || d1 < epsilon || d2 < epsilon || d3 < epsilon);
	//return (u >= 0) && (v >= 0) && (u + v <= 1);
}
bool Triangle::IsPointInTriangleRange(Vec3d  p, Vec3d a,Vec3d b,Vec3d c, float s)
{

	//-------------- scaling ------------Triangle::----
    /*
     *  UT_Vector3 v0 = c - a;
        UT_Vector3 v1 = b - a;
        UT_Vector3 v2 = p - a;

        // Compute dot products
        float dot00 = dot(v0, v0);
        float dot01 = dot(v0, v1);
        float dot02 = dot(v0, v2);
        float dot11 = dot(v1, v1);
        float dot12 = dot(v1, v2);

        // Compute barycentric coordinates
        float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
        float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

        // Check if point is in triangle
        return (u >= -epsilon) && (v >= -epsilon) && (u + v <= 1+epsilon);
     */
	
	
    float epsilon = s;
    Vec3d v0 = c - a;
    Vec3d v1 = b - a;
    Vec3d v2 = p - a;


    // Compute dot products
    float dot00 = v0.dot(v0);
    float dot01 = v0.dot(v1);
    float dot02 = v0.dot(v2);
    float dot11 = v1.dot(v1);
    float dot12 = v1.dot(v2);

    // Compute barycentric coordinates
    float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Check if point is in triangle
    return (u >= -epsilon) && (v >= -epsilon) && (u + v <= 1+epsilon);

	
}

Vec3d Triangle::GetBarycentricPosition(Vec3d A,Vec3d B, Vec3d C, Vec3d a, Vec3d b, Vec3d c, Vec3d position)
{
    /*
    double totalAir = TriangleAir(A,B,C);
    double apbAir = TriangleAir(A,B,position);
    double bpcAir = TriangleAir(B,C,position);
    double apcAir = TriangleAir(A,C, position);

    double alpha = bpcAir/totalAir;
    double beta =  apcAir/totalAir;
    double gamma = apbAir/totalAir;

	a *= alpha;
	b *= beta;
	c *= gamma;

	return a+b+c;
    */

    double a1 = ((B.y-C.y)*(position.x-C.x)+(C.x-B.x)*(position.y-C.y));
    double a2 = ((B.y-C.y)*(A.x-C.x)+(C.x-B.x)*(A.y-C.y));
    double b1 = ((C.y-A.y)*(position.x-C.x)+(A.x-C.x)*(position.y-C.y));
    double b2 = ((B.y-C.y)*(A.x-C.x)+(C.x-B.x)*(A.y-C.y));

    double alpha = a1/a2 ;
    double beta = b1/b2 ;
    double gamma = 1-alpha-beta;
    a *= alpha;
    b *= beta;
    c *= gamma;

    return a+b+c;


}




Vec3d Triangle::ProjectPointOnTriangle(Vec3d p,  Vec3d a,Vec3d b,Vec3d c)
{
	/*

	UT_Vector3 v1 = a - b;
	UT_Vector3 v2 = a - c;
	UT_Vector3 n = cross(v1,v2);GRID
    n.normalize();

	UT_Vector3 r0 = a;
	UT_Vector3 v3 = r0 - p;
	UT_Vector3 v3m = v3-n*dot(n,v3);
	UT_Vector3 rm = r0-v3m;
	return rm;
	*/
    Vec3d v1 = b - a;
    Vec3d v2 = c - a;
    Vec3d n = v1.cross(v2);
    n.normalize();

    Vec3d r0 = a;
    Vec3d v3 = p - r0;
    float dot = n.dot(v3);
    Vec3d ndot = n;
	ndot *= dot;
    Vec3d v3m = v3-ndot;
    Vec3d rm = r0 + v3m;
	return rm;
}


Vec3d Triangle::GetUVBarycenter(Vec3d a,Vec3d b, Vec3d c,Vec3d A,Vec3d B, Vec3d C, Vec3d P)
{
	/*    
	UT_Vector3 v0 = c - a;
	UT_Vector3 v1 = b - a;
	UT_Vector3 v2 = P - a;

	// Compute dot products
	float dot00 = dot(v0, v0);
	float dot01 = dot(v0, v1);
	float dot02 = dot(v0, v2);
	float dot11 = dot(v1, v1);
	float dot12 = dot(v1, v2);

	// Compute barycentric coordinatesTriangle::
	float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	return A + u*(C-A) + v* (B-A);
	*/
	
	// Compute vectors        
    //Vec3d v0 = c - a;
    //Vec3d v1 = b - a;
    //Vec3d v2 = P - a;

    Vec3d v0 = sous(a,c);
    Vec3d v1 = sous(a,b);
    Vec3d v2 = sous(a,P);


	
	// Compute dot products
    float dot00 = v0.dot(v0);
    float dot01 = v0.dot(v1);
    float dot02 = v0.dot(v2);
    float dot11 = v1.dot(v1);
    float dot12 = v1.dot(v2);

	// Compute barycentric coordinates
	float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	//return A + u*(C-A) + v* (B-A);
    Vec3d ac = sous(A,C);
    Vec3d ab = sous(A,B);
	ac *= u;
	ab *= v;
    Vec3d r;
	return A + ac + ab;
}

void Triangle::BoundingBox(Vec3d a,Vec3d b,Vec3d c, Vec3d &min,Vec3d &max)
{
	min.x = a.x;
	min.y = a.y;
	min.z = a.z;
	max = min;

	// check min

	if (min.x > b.x)
	{
		min.x = b.x;
	}
	if (min.y > b.y)
	{
		min.y = b.y;
	}
	if (min.z > b.z)
	{
		min.z = b.z;
	}

	if (min.x > c.x)
	{
		min.x = c.x;
	}
	if (min.y > c.y)
	{
		min.y = c.y;
	}
	if (min.z > c.z)
	{
		min.z = c.z;
	}

	//check max

	if (max.x < b.x)
	{
		max.x = b.x;
	}
	if (max.y < b.y)
	{
		max.y = b.y;
	}
	if (max.z < b.z)
	{
		max.z = b.z;
	}

	if (max.x < c.x)
	{
		max.x = c.x;
	}
	if (max.y < c.y)
	{
		max.y = c.y;
	}
	if (max.z < c.z)
	{
		max.z = c.z;
	}
}

void Triangle::BoundingBox2D(Vec3d a,Vec3d b,Vec3d c, Vec3d &min,Vec3d &max)
{
    min.x = a.x;
    min.y = a.y;
    //min.z = a.z;
    max = min;

    // check min

    if (min.x > b.x)
    {
        min.x = b.x;
    }
    if (min.y > b.y)
    {
        min.y = b.y;
    }

    if (min.x > c.x)
    {
        min.x = c.x;
    }
    if (min.y > c.y)
    {
        min.y = c.y;
    }

    //check max

    if (max.x < b.x)
    {
        max.x = b.x;
    }
    if (max.y < b.y)
    {
        max.y = b.y;
    }

    if (max.x < c.x)
    {
        max.x = c.x;
    }
    if (max.y < c.y)
    {
        max.y = c.y;
    }

}

