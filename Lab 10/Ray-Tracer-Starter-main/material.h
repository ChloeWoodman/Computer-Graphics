#pragma once

#include "common.h"
#include "geometry.h"
#include "hittable.h"

// Forward declaration of hit_record
struct hit_record;

// Abstract base class for materials
class material {
public:
    // Pure virtual function for computing scattered ray
    virtual bool scatter(const Ray& r_in, const hit_record& rec, Colour& attenuation, Ray& scattered) const = 0;
};

// Lambertian material (diffuse reflection)
class lambertian : public material {
public:
    // Constructor
    lambertian(const Colour& a) : albedo(a) {}

    // Implementation of scatter function for Lambertian material
    virtual bool scatter(
        const Ray& r_in, const hit_record& rec, Colour& attenuation, Ray& scattered
    ) const override {
        auto scatter_direction = rec.normal + Vec3f().random_in_unit_sphere();
        // Catch degenerate scatter direction
        if (scatter_direction.near_zero())
            scatter_direction = rec.normal;
        scattered = Ray(rec.p, scatter_direction);
        attenuation = albedo;
        return true;
    }

public:
    //albedo (colour) of material
    Colour albedo;
};

// Compute reflection direction for a given incident vector and surface normal
Vec3f reflect(const Vec3f& v, const Vec3f& n) {
    return v - 2 * v.dotProduct(n) * n;
}

// Metal material (specular reflection)
class metal : public material {
public:
    metal(const Colour& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}

    // Implementation of scatter function for metal material
    virtual bool scatter(
        const Ray& r_in, const hit_record& rec, Colour& attenuation, Ray& scattered
    ) const override {
        Vec3f reflected = reflect(r_in.direction().normalize(), rec.normal);
        scattered = Ray(rec.p, reflected + fuzz * Vec3f().random_in_unit_sphere());
        attenuation = albedo;
        return (scattered.direction().dotProduct(rec.normal) > 0);
    }

public:
    Colour albedo;
    double fuzz;
};
