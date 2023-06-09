#pragma once

#include "hittable.h"
#include "geometry.h"
#include "aabb.h"

// Declaration of the 'triangle' class, which is derived from the 'hittable' class
class triangle : public hittable {
public:
    // Default constructor for the 'triangle' class
    triangle() {}
    // Constructor for the 'triangle' class that takes in three vertices, and a shared pointer to a material
    triangle(Point3f vert0, Point3f vert1, Point3f vert2, shared_ptr<material> m) : v0(vert0), v1(vert1), v2(vert2), mat_ptr(m) {
        // Calculates the normal of the triangle using the cross product of two vectors
        normal = (v1 - v0).crossProduct(v2 - v0);
    };

    // Constructor for the 'triangle' class that takes in three vertices, a normal vector, and a shared pointer to a material
    triangle(Point3f vert0, Point3f vert1, Point3f vert2, Vec3f vn, shared_ptr<material> m) : v0(vert0), v1(vert1), v2(vert2), normal(vn), mat_ptr(m) {};

    triangle(Point3f vert0, Point3f vert1, Point3f vert2, Vec3f _vn0, Vec3f _vn1, Vec3f _vn2, shared_ptr<material> m) : v0(vert0), v1(vert1), v2(vert2), vn0(_vn0), vn1(_vn1), vn2(_vn2), mat_ptr(m) 
    {
        normal = vn0 + vn1 + vn2 / 3;
    };
    
    // Overrides the 'hit' function of the 'hittable' class
    virtual bool hit(
        const Ray& r, double t_min, double t_max, hit_record& rec) const override;

    // Overrides the 'bounding_box' function of the 'hittable' class
    virtual bool bounding_box(aabb& output_box) const override;
private:
    // Utility function to map points on a unit triangle to UV coordinates
    static void get_tri_uv(const Point3f& p, double& u, double& v) {
        // p: a given point on the unit triangle.
        // u: returned value [0,1] of the x-coordinate of p.
        // v: returned value [0,1] of the y-coordinate of p.

        u = p.x;
        v = p.y;
    }

public:
    // The three vertices of the triangle
    Point3f v0, v1, v2;
    // The normal vector of the triangle
    Vec3f normal;
    Vec3f vn0, vn1, vn2;
    // A shared pointer to the material of the triangle
    shared_ptr<material> mat_ptr;
};

//looking just for min and max values for the x, y and z per vertex of triangle
//...these define the box
//bounding box must be implemented for each primitive due to deriving from hittable class
inline bool triangle::bounding_box(aabb& output_box) const
{
    float min[3];
    float max[3];
    for (int i = 0; i < 3; i++) { // For each dimension
    // Calculates the minimum and maximum values of the vertices in the triangle
        min[i] = std::min(v0[i], std::min(v1[i], v2[i]));
        max[i] = std::max(v0[i], std::max(v1[i], v2[i]));
    }
    output_box = aabb(Vec3f(min[0], min[1], min[2]), Vec3f(max[0], max[1], max[2])); // Sets the output box to the calculated AABB
    return true;
}


// ray-triangle intersection: 
// https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution
// be aware, faster methods exist: 
// https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection

bool triangle::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
    double thit, t, u, v;

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

    if (t < t_min || t > t_max) return false;

    rec.t = t;
    rec.p = r.at(t);
    //rec.normal = normal;
    if (vn0.x > 1e-3) {
        rec.normal = vn1 * u + vn2 * v + (1.0f - u - v) * vn0;
    }
    else {rec.normal = normal;}
    rec.mat_ptr = mat_ptr;

    // Calculate and set UV coordinates
    Point3f p_in_triangle = (rec.p - v0) / v0v1.length();
    get_tri_uv(p_in_triangle, rec.u, rec.v);

    // Calculate and set normal
    //rec.normal = normal;

    return true;
}