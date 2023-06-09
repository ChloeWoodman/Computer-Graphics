// This code defines a camera class that allows us to generate rays through pixels in the image plane
#pragma once

#include "common.h"

class camera {
public:
    // Constructor sets up the camera's parameters based on the desired aspect ratio and viewport height
    camera() {
        auto aspect_ratio = 16.0 / 9.0;
        auto viewport_height = 2.0;
        auto viewport_width = aspect_ratio * viewport_height;
        auto focal_length = 1.0;

        origin = Point3f(0, 0, 0);
        horizontal = Vec3f(viewport_width, 0.0, 0.0);
        vertical = Vec3f(0.0, viewport_height, 0.0);
        lower_left_corner = origin - horizontal / 2 - vertical / 2 - Vec3f(0, 0, focal_length);
    }

    // Generates a ray through a given pixel (u,v) in the image plane
    Ray get_ray(double u, double v) const {
        return Ray(origin, lower_left_corner + u * horizontal + v * vertical - origin);
    }

private:
    Point3f origin; // Camera position
    Point3f lower_left_corner; // Lower left corner of the image plane
    Vec3f horizontal; // Horizontal vector of the image plane
    Vec3f vertical; // Vertical vector of the image plane
};