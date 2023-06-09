#pragma once
#include "ray.h"
#include "common.h"
#include "aabb.h"

// forward declaration of the material class
class material;

// a struct to store information about the intersection
struct hit_record {
    Point3f p; // intersection point
    Vec3f normal; // normal at the intersection point
    double t; // parameter value for the ray at the intersection point
    bool front_face; // true if the ray hit the front face of the object

    // a utility function to set the front face and normal
    inline void set_face_normal(const Ray& r, const Vec3f& outward_normal) {
        front_face = (r.direction().dotProduct(outward_normal)) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }

    shared_ptr<material> mat_ptr; // pointer to the material of the object
};

// an abstract class to define a hittable object
class hittable {
public:
    // pure virtual, must be overwritten in derived classes
    virtual bool hit(const Ray& r, double t_min, double t_max, hit_record& rec) const = 0;
    virtual bool bounding_box(aabb& output_box) const = 0; //compute bounding boxes of all hittables to comprise the scene
};
