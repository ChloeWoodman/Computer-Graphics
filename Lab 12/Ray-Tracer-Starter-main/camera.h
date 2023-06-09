// This code defines a camera class that allows us to generate rays through pixels in the image plane
#pragma once

#include "common.h"

class camera {
public:
    // Constructor sets up the camera's parameters based on the desired aspect ratio and viewport height
    camera(
        Point3f lookfrom,
        Point3f lookat,
        Vec3f   vup,
        double vfov, // vertical field-of-view in degrees
        double aspect_ratio,
        double aperture,
        double focus_dist
    ) {
        auto theta = degrees_to_radians(vfov);
        auto h = tan(theta / 2);
        auto viewport_height = 2.0 * h;
        auto viewport_width = aspect_ratio * viewport_height;
        w = (lookfrom - lookat).normalize();
        u = (vup.crossProduct(w)).normalize();
        v = w.crossProduct(u);
        origin = lookfrom;
        horizontal = focus_dist * viewport_width * u;
        vertical = focus_dist * viewport_height * v;
        lower_left_corner = origin - horizontal / 2 - vertical / 2 - focus_dist * w;
        lens_radius = aperture / 2;
    }
    Ray get_ray(double s, double t) const {
        Vec3f rd = lens_radius * Vec3f().random_in_unit_disk();
        Vec3f offset = u * rd.x + v * rd.y;
        return Ray(
            origin + offset,
            lower_left_corner + s * horizontal + t * vertical - origin - offset
        );
    }

private:
    Point3f origin; // Camera position
    Point3f lower_left_corner; // Lower left corner of the image plane
    Vec3f horizontal; // Horizontal vector of the image plane
    Vec3f vertical; // Vertical vector of the image plane
    Vec3f u, v, w;
    double lens_radius;
};