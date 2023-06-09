#pragma once

#include "hittable.h"

#include <memory>
#include <vector>

using std::shared_ptr; //using shared pointer from memory management library
using std::make_shared; //using make_shared from memory management library

class hittable_list : public hittable { //class declaration for hittable_list, which inherits from hittable
public:
    // Default constructor
    hittable_list() {} //default constructor

    // Constructor that adds an object to the list
    hittable_list(shared_ptr<hittable> object) { add(object); }

    // Clear the list of objects
    void clear() { objects.clear(); } 

    // Add an object to the list
    void add(shared_ptr<hittable> object) { objects.push_back(object); } //function to add object to list

    // Check if a ray hits any objects in the list and return the closest hit
    virtual bool hit(const Ray& r, double t_min, double t_max, hit_record& rec) const override; 

    //function to calculate the bounding box of the hittable_list
    virtual bool bounding_box(aabb& output_box) const override;

public:
    // A vector to store pointers to the hittable objects in the list
    std::vector<shared_ptr<hittable>> objects;
};

//since this class lists object primitives, bounding box can be stored at contrustion of the list
inline bool hittable_list::bounding_box(aabb& output_box) const
{
    // If there are no objects in the list, return false
    if (objects.empty()) return false;

    // Initialize temporary variables
    aabb temp_box;
    bool first_box = true;

    // Iterate through all objects in the list
    for (const auto& object : objects) {
        // If an object does not have a bounding box, return false
        if (!object->bounding_box(temp_box)) return false;

        // Update the output_box to include the bounding box of the current object
        output_box = first_box ? temp_box : surrounding_box(output_box, temp_box);
        first_box = false;
    }

    // Return true to indicate that all objects have a bounding box
    return true;
}


// Implementation of the hit function for hittable_list class
// This function checks if a ray hits any of the objects in the list and updates
// the hit_record parameter passed as reference if it does
bool hittable_list::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
    // Initialize temporary variables
    hit_record temp_rec;
    bool hit_anything = false;
    auto closest_so_far = t_max;

    // Iterate through all objects in the list
    for (const auto& object : objects) {
        // Check if the ray hits the object within the given time range and update
        // the hit_record if it does
        if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }

    // Return true if the ray hits any object in the list
    return hit_anything;
}

class translate : public hittable {
public:
    translate(shared_ptr<hittable> p, const Vec3f& displacement)
        : ptr(p), offset(displacement) {}

    virtual bool hit(
        const Ray& r, double t_min, double t_max, hit_record& rec) const override;

    virtual bool bounding_box(aabb& output_box) const override;

public:
    shared_ptr<hittable> ptr;
    Vec3f offset;
};

bool translate::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
    Ray moved_r(r.origin() - offset, r.direction(), r.time());
    if (!ptr->hit(moved_r, t_min, t_max, rec))
        return false;

    rec.p.x + offset.x;
    rec.p.y + offset.y;
    rec.p.z + offset.z;
    rec.set_face_normal(moved_r, rec.normal);

    return true;
}

bool translate::bounding_box(aabb& output_box) const {
    if (!ptr->bounding_box(output_box))
        return false;

    output_box = aabb(
        output_box.min() + offset,
        output_box.max() + offset);

    return true;
}

class rotate_y : public hittable {
public:
    rotate_y(shared_ptr<hittable> p, double angle);

    virtual bool hit(
        const Ray& r, double t_min, double t_max, hit_record& rec) const override;

    virtual bool bounding_box(aabb& output_box) const override {
        output_box = bbox;
        return hasbox;
    }

public:
    shared_ptr<hittable> ptr;
    double sin_theta;
    double cos_theta;
    bool hasbox;
    aabb bbox;
};

rotate_y::rotate_y(shared_ptr<hittable> p, double angle) : ptr(p) {
    auto radians = degrees_to_radians(angle);
    sin_theta = sin(radians);
    cos_theta = cos(radians);
    hasbox = ptr->bounding_box(bbox);

    Point3f min(infinity, infinity, infinity);
    Point3f max(-infinity, -infinity, -infinity);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                auto x = i * bbox.max().x + (1 - i) * bbox.min().x;
                auto y = j * bbox.max().y + (1 - j) * bbox.min().y;
                auto z = k * bbox.max().z + (1 - k) * bbox.min().z;

                auto newx = cos_theta * x + sin_theta * z;
                auto newz = -sin_theta * x + cos_theta * z;

                Vec3f tester(newx, y, newz);

                for (int c = 0; c < 3; c++) {
                    min[c] = fmin(min[c], tester[c]);
                    max[c] = fmax(max[c], tester[c]);
                }
            }
        }
    }

    bbox = aabb(min, max);
}

bool rotate_y::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
    auto origin = r.origin();
    auto direction = r.direction();

    origin[0] = cos_theta * r.origin()[0] - sin_theta * r.origin()[2];
    origin[2] = sin_theta * r.origin()[0] + cos_theta * r.origin()[2];

    direction[0] = cos_theta * r.direction()[0] - sin_theta * r.direction()[2];
    direction[2] = sin_theta * r.direction()[0] + cos_theta * r.direction()[2];

    Ray rotated_r(origin, direction, r.time());

    if (!ptr->hit(rotated_r, t_min, t_max, rec))
        return false;

    auto p = rec.p;
    auto normal = rec.normal;

    p[0] = cos_theta * rec.p[0] + sin_theta * rec.p[2];
    p[2] = -sin_theta * rec.p[0] + cos_theta * rec.p[2];

    normal[0] = cos_theta * rec.normal[0] + sin_theta * rec.normal[2];
    normal[2] = -sin_theta * rec.normal[0] + cos_theta * rec.normal[2];

    rec.p = p;
    rec.set_face_normal(rotated_r, normal);

    return true;
}