#pragma once
#include "common.h"
#include "hittable.h"
#include "aabb.h"
#include <cmath>


// Definition of the moving_sphere class
class moving_sphere : public hittable {
public:
	moving_sphere() {}
	moving_sphere(Vec3f cen0, Vec3f cen1, double _time0, double _time1, double r, shared_ptr<material> m)
		: centre0(cen0), centre1(cen1), time0(_time0), time1(_time1), radius(r), mat_ptr(m) {};

	virtual bool hit(const Ray& r, double t_min, double t_max, hit_record& rec) const override;

	// Function for calculating the centre of the sphere at a given time
	Point3f centre(double time) const;

	// Function for computing a bounding box
	virtual bool bounding_box(aabb& output_box) const override;

public:
	Vec3f centre0, centre1;
	double time0, time1;
	double radius;
	shared_ptr<material> mat_ptr;
};

// Implementation of centre function
Point3f moving_sphere::centre(double time) const {
	return centre0 + ((time - time0) / (time1 - time0)) * (centre1 - centre0);
}

// Implementation of bounding_box function
inline bool moving_sphere::bounding_box(aabb& output_box) const {
	output_box = aabb(centre0 - Vec3f(radius, radius, radius), centre1 + Vec3f(radius, radius, radius));
	return true;
}

// Implementation of hit function
bool moving_sphere::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
	Vec3f oc = r.origin() - centre(r.time());
	auto a = r.direction().norm();
	auto half_b = oc.dotProduct(r.direction());
	auto c = oc.norm() - radius * radius;
	auto discriminant = half_b * half_b - a * c;
	if (discriminant < 0) return false;
	auto sqrtd = sqrt(discriminant);

	// Find the nearest root that lies in the acceptable range
	auto root = (-half_b - sqrtd) / a;
	if (root < t_min || t_max < root) {
		root = (-half_b + sqrtd) / a;
		if (root < t_min || t_max < root)
			return false;
	}

	// Update hit record data accordingly
	rec.t = root;
	rec.p = r.at(rec.t);
	auto outward_normal = (rec.p - centre(r.time())) / radius;
	rec.set_face_normal(r, outward_normal);
	rec.mat_ptr = mat_ptr;

	return true;
}