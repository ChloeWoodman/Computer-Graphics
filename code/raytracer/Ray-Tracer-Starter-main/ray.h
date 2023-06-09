#pragma once
#include "geometry.h"

class Ray {
public:
    Ray() {}
    Ray(const Point3f& origin, const Vec3f& direction, double time = 0.0)
        : orig(origin), dir(direction), tm(time)
    {}
    //Ray(const Point3f& origin, const Vec3f& direction, double time = 0.0)
    //    : orig(origin), dir(direction)
    //{}

    Point3f origin() const { return orig; }
    Vec3f direction() const { return dir; }
    double time() const { return tm; }

    Point3f at(double t) const {
        return orig + t * dir;
    }

public:
    Point3f orig;
    Vec3f dir;
    double tm;
};