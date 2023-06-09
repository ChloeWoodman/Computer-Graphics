#pragma once
#include "hittable.h"
#include "geometry.h"

// Sphere class inheriting from hittable class
class Sphere : public hittable {
public:
	// Constructors
	Sphere() {}
	Sphere(Vec3f cen, double r, shared_ptr<material> m) : centre(cen), radius(r), mat_ptr(m) {};

	// Override hit function from hittable class
	virtual bool hit(const Ray& r, double t_min, double t_max, hit_record& rec) const override;
	virtual bool bounding_box(aabb& output_box) const override;

public:
	// Sphere properties
	Vec3f centre;
	double radius;
	shared_ptr<material> mat_ptr;
};

//bounding box must be implemented for each primitive due to deriving from hittable class
inline bool Sphere::bounding_box(aabb& output_box) const
{
	output_box = aabb(centre - Vec3f(radius, radius, radius), centre + Vec3f(radius, radius, radius));
	return true;
}

// Definition of hit function
bool Sphere::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
	Vec3f oc = r.origin() - centre;
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
	Vec3f outward_normal = (rec.p - centre) / radius;
	rec.set_face_normal(r, outward_normal);
	rec.mat_ptr = mat_ptr;

	return true;
}