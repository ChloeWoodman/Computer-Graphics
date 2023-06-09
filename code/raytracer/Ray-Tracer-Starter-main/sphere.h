#pragma once
#include "hittable.h"
#include "geometry.h"
#include "aabb.h"
#include <cmath>

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

private:
	// Utility function to map points on a unit sphere to UV coordinates
	static void get_sphere_uv(const Point3f& p, double& u, double& v) {
		// p: a given point on the sphere of radius one, centered at the origin.
		// u: returned value [0,1] of angle around the Y axis from X=-1.
		// v: returned value [0,1] of angle from Y=-1 to Y=+1.
		//     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
		//     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
		//     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>

		auto theta = acos(-p.y);
		auto phi = atan2(-p.z, p.x) + pi;

		u = phi / (2 * pi);
		v = theta / pi;
	}
};

// Bounding box function for Sphere
inline bool Sphere::bounding_box(aabb& output_box) const
{
	// Calculate bounding box
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

	// Check if there is a hit
	if (discriminant < 0) return false;

	// Calculate the roots of the quadratic equation
	auto sqrtd = sqrt(discriminant);
	auto root = (-half_b - sqrtd) / a;

	// Check if the root is within the acceptable range
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
	get_sphere_uv(outward_normal, rec.u, rec.v);
	rec.mat_ptr = mat_ptr;

	return true;
}