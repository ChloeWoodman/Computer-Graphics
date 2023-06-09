#pragma once

#include "hittable.h"

#include <memory>
#include <vector>

using std::shared_ptr;
using std::make_shared;

class hittable_list : public hittable {
public:
    // Default constructor
    hittable_list() {}

    // Constructor that adds an object to the list
    hittable_list(shared_ptr<hittable> object) { add(object); }

    // Clear the list of objects
    void clear() { objects.clear(); }

    // Add an object to the list
    void add(shared_ptr<hittable> object) { objects.push_back(object); }

    // Check if a ray hits any objects in the list and return the closest hit
    virtual bool hit(const Ray& r, double t_min, double t_max, hit_record& rec) const override;

public:
    // A vector to store pointers to the hittable objects in the list
    std::vector<shared_ptr<hittable>> objects;
};

// Implementation of the hit function for hittable_list class
bool hittable_list::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
    hit_record temp_rec;
    bool hit_anything = false;
    auto closest_so_far = t_max;

    // Iterate through all objects in the list
    for (const auto& object : objects) {
        // Check if the ray hits the object and update the hit record if it does
        if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }

    return hit_anything;

}