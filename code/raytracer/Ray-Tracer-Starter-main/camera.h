#pragma once

#include "common.h"

class camera {
public:
    camera(
        Point3f lookfrom,
        Point3f lookat,
        Vec3f vup,
        double vfov,
        double aspect_ratio,
        double aperture,
        double focus_dist,
        double _time0 = 0,
        double _time1 = 0
    ) {
        updateCamera(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, focus_dist, _time0, _time1);
    }

    void updateCamera(
        Point3f lookfrom,
        Point3f lookat,
        Vec3f vup,
        double vfov,
        double aspect_ratio,
        double aperture,
        double focus_dist,
        double _time0,
        double _time1
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
        time0 = _time0;
        time1 = _time1;
    }

    Ray get_ray(double s, double t) const {
        Vec3f rd = lens_radius * Vec3f().random_in_unit_disk();
        Vec3f offset = u * rd.x + v * rd.y;
        return Ray(
            origin + offset,
            lower_left_corner + s * horizontal + t * vertical - origin - offset,
            random_double(time0, time1)
        );
    }

private:
    Point3f origin;
    Point3f lower_left_corner;
    Vec3f horizontal;
    Vec3f vertical;
    Vec3f u, v, w;
    double lens_radius;
    double time0, time1;
};
