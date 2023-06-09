#pragma once
#include "common.h"
#include "hittable.h"
#include <cmath>
#include "aabb.h"


// Define the class for a moving triangle that inherits from the hittable class
class moving_triangle : public hittable {
public:
	// Default constructor
	moving_triangle() {}
	// Constructor that initializes the triangle's vertices, times, and material pointer
	moving_triangle(
		Point3f vert0, Point3f vert1, Point3f vert2, double _time0, double _time1, shared_ptr<material> m)
		: v0(vert0), v1(vert1), v2(vert2), time0(_time0), time1(_time1), mat_ptr(m) {
		normal = (v1 - v0).crossProduct(v2 - v0); // Calculate the triangle's normal
	};

	// Constructor that initializes the triangle's vertices, times, normal, and material pointer
	moving_triangle(Point3f vert0, Point3f vert1, Point3f vert2, double _time0, double _time1, Vec3f vn, shared_ptr<material> m)
		: v0(vert0), v1(vert1), v2(vert2), time0(_time0), time1(_time1), normal(vn), mat_ptr(m) {};

	// Override the hit function of the hittable class
	virtual bool hit(const Ray& r, double t_min, double t_max, hit_record& rec) const override;

	// Calculate the centre of the triangle at a given time
	Point3f centre(double time) const;

	// Override the bounding_box function of the hittable class to define the bounding box of the triangle
	virtual bool bounding_box(aabb& output_box) const override;

public:
	Vec3f vert0, vert1;
	double time0, time1;
	Point3f v0, v1, v2;
	Vec3f normal;
	shared_ptr<material> mat_ptr;
};

// Calculate the centre of the triangle at a given time
Point3f moving_triangle::centre(double time) const {
	return vert0 + ((time - time0) / (time1 - time0)) * (vert1 - vert0);
}

// Override the bounding_box function of the hittable class to define the bounding box of the triangle
inline bool moving_triangle::bounding_box(aabb& output_box) const
{
	float min[3];
	float max[3];
	for (int i = 0; i < 3; i++) { // for each dimension
	// calculate minimum and maximum values of the vertices in the triangle
		min[i] = std::min(v0[i], std::min(v1[i], v2[i]));
		max[i] = std::max(v0[i], std::max(v1[i], v2[i]));
	}
	output_box = aabb(Vec3f(min[0], min[1], min[2]), Vec3f(max[0], max[1], max[2])); // Define the bounding box
	return true;
}

// Definition of hit function
bool moving_triangle::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {

	// Variables for storing intersection information
	float thit, t, u, v;

	// Calculate vectors for edges and the normal of the triangle
	Vec3<float> v0v1 = v1 - v0;
	Vec3<float> v0v2 = v2 - v0;

	// Calculate the cross product of the ray direction and v0v2
	Vec3<float> pvec = r.direction().crossProduct(v0v2);

	// Calculate the determinant of the matrix
	float det = pvec.dotProduct(v0v1);
	float kEpsilon = 0.00001;

	// If the determinant is negative, the triangle is backfacing
	// If the determinant is close to 0, the ray misses the triangle
	if (det < kEpsilon) return false;

	// Calculate inverse of the determinant
	float invDet = 1 / det;

	// Calculate u parameter and test bounds
	Vec3<float> tvec = r.origin() - v0;
	u = tvec.dotProduct(pvec) * invDet;
	if (u < 0 || u > 1) return false;

	// Calculate v parameter and test bounds
	Vec3<float> qvec = tvec.crossProduct(v0v1);
	v = r.direction().dotProduct(qvec) * invDet;
	if (v < 0 || u + v > 1) return false;

	// Calculate t parameter
	t = v0v2.dotProduct(qvec) * invDet;

	// If t is less than 0, the intersection is behind the camera
	if (t < 0) return false;

	// Set hit record data
	rec.p = r.at(t);
	rec.t = t;
	// TODO: fix normal calculation as materials depend on this. How do we determine correct triangle normals?
	rec.normal = normal;
	rec.mat_ptr = mat_ptr;

	return true;
}
