#pragma once
#include "common.h"
#include "hittable.h"
#include "material.h"
#include "texture.h"

// constant_medium is a class that inherits from the hittable class. It represents a volume with constant density that a ray can hit.
class constant_medium : public hittable {
public:
    // Constructor: takes a shared_ptr to a hittable object, a density, and a shared_ptr to a texture object
    constant_medium(shared_ptr<hittable> b, double d, shared_ptr<texture> a)
        : boundary(b),
        neg_inv_density(-1 / d), // calculate the negative inverse density
        phase_function(make_shared<isotropic>(a))
    {}

    // Constructor: takes a shared_ptr to a hittable object, a density, and a Colour object
    constant_medium(shared_ptr<hittable> b, double d, Colour c)
        : boundary(b),
        neg_inv_density(-1 / d), // calculate the negative inverse density
        phase_function(make_shared<isotropic>(c))
    {}

    // Overriding hit method from hittable class
    virtual bool hit(
        const Ray& r, double t_min, double t_max, hit_record& rec) const override;

    // Overriding bounding_box method from hittable class
    virtual bool bounding_box(aabb& output_box) const override {
        return boundary->bounding_box(output_box);
    }

public:
    // Public member variables
    shared_ptr<hittable> boundary; // a boundary that can be hit
    shared_ptr<material> phase_function; // the material of the medium
    double neg_inv_density; // the negative inverse density of the medium
};

// Implementation of the hit method for the constant_medium class
bool constant_medium::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
    // Debugging flags
    const bool enableDebug = false;
    const bool debugging = enableDebug && random_double() < 0.00001;

    hit_record rec1, rec2; // two hit records

    // Check if the ray hits the boundary
    if (!boundary->hit(r, -infinity, infinity, rec1))
        return false;

    if (!boundary->hit(r, rec1.t + 0.0001, infinity, rec2))
        return false;

    if (debugging) std::cerr << "\nt_min=" << rec1.t << ", t_max=" << rec2.t << '\n';

    // Adjust the hit points based on t_min and t_max
    if (rec1.t < t_min) rec1.t = t_min;
    if (rec2.t > t_max) rec2.t = t_max;

    if (rec1.t >= rec2.t)
        return false;

    if (rec1.t < 0)
        rec1.t = 0;

    // Calculate the length of the ray, the distance inside the boundary, and the hit distance
    const auto ray_length = r.direction().length();
    const auto distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
    const auto hit_distance = neg_inv_density * log(random_double());

    // If the hit_distance is greater than the distance inside the boundary, return false
    if (hit_distance > distance_inside_boundary)
        return false;

    // Set the hit record's t value and p value
    rec.t = rec1.t + hit_distance / ray_length;
    rec.p = r.at(rec.t);

    if (debugging) {
        std::cerr << "hit_distance = " << hit_distance << '\n'
            << "rec.t = " << rec.t << '\n'
            << "rec.p = " << rec.p << '\n';
    }

    // Set arbitrary normal vector, front_face, and material pointer
    rec.normal = Vec3f(1, 0, 0);  // arbitrary
    rec.front_face = true;     // also arbitrary
    rec.mat_ptr = phase_function;

    return true;
}
