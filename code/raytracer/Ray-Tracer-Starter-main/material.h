#pragma once

#include "common.h"
#include "geometry.h"
#include "hittable.h"
#include "texture.h"

// Forward declaration of hit_record structure
struct hit_record;

// Abstract base class for materials
class material {
public:
    // Pure virtual function for computing scattered ray
    virtual bool scatter(const Ray& r_in, const hit_record& rec, Colour& attenuation, Ray& scattered) const = 0;
    virtual Colour emitted() const { // Virtual function to compute emitted colour
        return Colour(0, 0, 0); // Default implementation returns black
    }
    virtual Colour emitted(double u, double v, const Point3f& p) const {
        return Colour(0, 0, 0);
    }
};

// Lambertian material (diffuse reflection)
class lambertian : public material {
public:
    lambertian(const Colour& a) : albedo(make_shared<solid_color>(a)) {}
    lambertian(shared_ptr<texture> a) : albedo(a) {}
    // Implementation of scatter function for Lambertian material
    virtual bool scatter(
        const Ray& r_in, const hit_record& rec, Colour& attenuation, Ray& scattered
    ) const override {
        auto scatter_direction = rec.normal + Vec3f().random_in_unit_sphere(); // Compute random scatter direction
        // Catch degenerate scatter direction
        if (scatter_direction.near_zero())
            scatter_direction = rec.normal;
        scattered = Ray(rec.p, scatter_direction, r_in.time()); // Update scattered ray
        //scattered = Ray(rec.p, scatter_direction);
        attenuation = albedo->value(rec.u, rec.v, rec.p);
        //attenuation = albedo; // Update attenuation
        return true; // Return true to indicate scattering occurred
    }
public:
    shared_ptr<texture> albedo;
    //Colour albedo; // Albedo (colour) of the material
};

// Compute reflection direction for a given incident vector and surface normal
Vec3f reflect(const Vec3f& v, const Vec3f& n) {
    return v - 2 * v.dotProduct(n) * n; // Compute and return the reflected vector
}

// Metal material (specular reflection)
class metal : public material {
public:
    metal(const Colour& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}
    // Implementation of scatter function for metal material
    virtual bool scatter(
        const Ray& r_in, const hit_record& rec, Colour& attenuation, Ray& scattered
    ) const override {
        Vec3f reflected = reflect(r_in.direction().normalize(), rec.normal); // Compute reflected vector
        scattered = Ray(rec.p, reflected + fuzz * Vec3f().random_in_unit_sphere(), r_in.time()); // Update scattered ray
        //scattered = Ray(rec.p, reflected + fuzz * Vec3f().random_in_unit_sphere());
        attenuation = albedo; // Update attenuation
        return (scattered.direction().dotProduct(rec.normal) > 0); // Return true if scattered ray is in the same hemisphere as the normal
    }
public:
    Colour albedo; // Albedo (colour) of the material
    double fuzz; // Fuzziness of the reflection
};

// Compute the refracted ray direction for a given incident ray, surface normal, and ratio of indices of refraction
Vec3f refract(const Vec3f& uv, const Vec3f& n, double etai_over_etat) {
    auto cos_theta = fmin(-uv.dotProduct(n), 1.0);
    Vec3f r_out_perp = etai_over_etat * (uv + cos_theta * n);
    Vec3f r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.norm())) * n;
    return r_out_perp + r_out_parallel;
}

// A dielectric material that refracts or reflects rays depending on the angle of incidence
class dielectric : public material {
public:
    // Constructor taking the index of refraction as an argument
    dielectric(double index_of_refraction) : ir(index_of_refraction) {}

    // Scatter function that returns true if a ray is scattered, false otherwise
    virtual bool scatter(
        const Ray& r_in, const hit_record& rec, Colour& attenuation, Ray& scattered
    ) const override {
        // Set the attenuation to white
        attenuation = Colour(1.0, 1.0, 1.0);

        // Calculate the ratio of indices of refraction
        double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;

        // Calculate the unit direction of the incident ray
        Vec3f unit_direction = r_in.direction().normalize();

        // Calculate the refracted ray direction
        Vec3f refracted = refract(unit_direction, rec.normal, refraction_ratio);

        // Calculate the cosine of the angle between the incident ray direction and the surface normal
        double cos_theta = fmin(-unit_direction.dotProduct(rec.normal), 1.0);

        // Calculate the sine of the angle between the incident ray direction and the surface normal
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

        // Check if the ray cannot be refracted
        bool cannot_refract = refraction_ratio * sin_theta > 1.0;

        // Declare a variable for the scattered ray direction
        Vec3f direction;

        // If the ray cannot be refracted or the reflectance is greater than a random number, reflect the ray
        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
            direction = reflect(unit_direction, rec.normal);
        // Otherwise, refract the ray
        else
            direction = refract(unit_direction, rec.normal, refraction_ratio);

        // Set the scattered ray to the reflected or refracted direction
        scattered = Ray(rec.p, direction, r_in.time());
        //scattered = Ray(rec.p, direction);

        // Return true, indicating that the ray has been scattered
        return true;
    }

public:
    // Index of refraction of the material
    double ir;
private:
    // Function to compute the reflectance given the cosine of the angle and the index of refraction
    static double reflectance(double cosine, double ref_idx) {
        // Use Schlick's approximation for reflectance.
        auto r0 = (1 - ref_idx) / (1 + ref_idx);
        r0 = r0 * r0;
        return r0 + (1 - r0) * pow((1 - cosine), 5);
    }
};

// A material that represents a diffuse light source
class diffuse_light : public material {
public:
    // Constructor accepting a texture that represents the light emission
    diffuse_light(shared_ptr<texture> a) : emit(a) {}

    // Constructor accepting a color that represents the light emission
    diffuse_light(Colour c) : emit(make_shared<solid_color>(c)) {}

    // This method checks if a ray scatters after hitting the material.
    // For a light source, this always returns false as light does not scatter further.
    virtual bool scatter(
        const Ray& r_in, const hit_record& rec, Colour& attenuation, Ray& scattered
    ) const override {
        return false;
    }

    // This method returns the color of the light that is being emitted.
    virtual Colour emitted(double u, double v, const Point3f& p) const override {
        return emit->value(u, v, p);
    }

public:
    // Emissive texture
    shared_ptr<texture> emit;
};


//class diffuse_light : public material {
//public:
//    diffuse_light() {}
//    diffuse_light(Colour c) : emit(make_shared<Colour>(c)) {}
//    virtual bool scatter(const Ray& r_in, const hit_record& rec, Colour& attenuation, Ray& scattered) const override {
//        return false;
//    }
//    virtual Colour emitted() const override {
//        return *emit;
//    }
//public:
//    shared_ptr<Colour> emit;
//};

// This is an isotropic material, meaning it scatters light equally in all directions.
class isotropic : public material {
public:
    // Constructor accepting a color that represents the material's albedo (reflectivity)
    isotropic(Colour c) : albedo(make_shared<solid_color>(c)) {}

    // Constructor accepting a texture that represents the material's albedo
    isotropic(shared_ptr<texture> a) : albedo(a) {}

    // This method checks if a ray scatters after hitting the material.
    // For an isotropic material, it always scatters in a random direction within a unit sphere.
    virtual bool scatter(
        const Ray& r_in, const hit_record& rec, Colour& attenuation, Ray& scattered
    ) const override {
        scattered = Ray(rec.p, Vec3f().random_in_unit_sphere(), r_in.time());
        attenuation = albedo->value(rec.u, rec.v, rec.p);
        return true;
    }

public:
    // Albedo texture
    shared_ptr<texture> albedo;
};