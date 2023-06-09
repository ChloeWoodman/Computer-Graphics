#pragma once

#include "hittable.h"
#include "geometry.h"

class triangle : public hittable {
public:
    triangle() {}
    triangle(Point3f vert0, Point3f vert1, Point3f vert2, shared_ptr<material> m) : v0(vert0), v1(vert1), v2(vert2), mat_ptr(m) {
        normal = (v1 - v0).crossProduct(v2 - v0);
    };
    triangle(Point3f vert0, Point3f vert1, Point3f vert2, Vec3f vn, shared_ptr<material> m) : v0(vert0), v1(vert1), v2(vert2), normal(vn), mat_ptr(m) {};

    virtual bool hit(
        const Ray& r, double t_min, double t_max, hit_record& rec) const override;
    virtual bool bounding_box(aabb& output_box) const override; //repeating what was done to sphere to here too
public:
    Point3f v0, v1, v2;
    Vec3f normal;
    shared_ptr<material> mat_ptr;
};

//looking just for min and max values for the x, y and z per vertex of triangle
//...these define the box
//bounding box must be implemented for each primitive due to deriving from hittable class
inline bool triangle::bounding_box(aabb& output_box) const
{
    float min[3];
    float max[3];
    for (int i = 0; i < 3; i++) { // for each dimension
    // calculate minimum and maximum values of the vertices in the triangle
        min[i] = std::min(v0[i], std::min(v1[i], v2[i]));
        max[i] = std::max(v0[i], std::max(v1[i], v2[i]));
    }
    output_box = aabb(Vec3f(min[0], min[1], min[2]), Vec3f(max[0], max[1], max[2]));
    return true;
}


// ray-triangle intersection: 
// https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution
// be aware, faster methods exist: 
// https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection
bool triangle::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {

    float thit, t, u, v;

    Vec3<float> v0v1 = v1 - v0;
    Vec3<float> v0v2 = v2 - v0;

    Vec3<float> pvec = r.direction().crossProduct(v0v2);

    float det = pvec.dotProduct(v0v1);
    float kEpsilon = 0.00001;

    // if the determinant is negative the triangle is backfacing
    // if the determinant is close to 0, the ray misses the triangle
    if (det < kEpsilon) return false;

    float invDet = 1 / det;

    Vec3<float> tvec = r.origin() - v0;
    u = tvec.dotProduct(pvec) * invDet;

    if (u < 0 || u > 1) return false;

    Vec3<float> qvec = tvec.crossProduct(v0v1);
    v = r.direction().dotProduct(qvec) * invDet;
    if (v < 0 || u + v > 1) return false;

    t = v0v2.dotProduct(qvec) * invDet;

    if (t < 0) return false; // intersection behind the camera

    rec.p = r.at(t);
    rec.t = t;
    // TODO: fix normal calculation as materials depend on this. How do we determine correct triangle normals?
    rec.normal = normal;
    rec.mat_ptr = mat_ptr;

    return true;
}

