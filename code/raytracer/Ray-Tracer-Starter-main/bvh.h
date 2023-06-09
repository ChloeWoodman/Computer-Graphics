#pragma once
#include "common.h"
#include "hittable.h"
#include "hittable_list.h"
#include <algorithm>
#include "aabb.h"

// Class declaration for BVH node which is a type of hittable object
class bvh_node : public hittable {
public:
    // Default constructor
    bvh_node();
    // Constructor that takes a list of hittable objects
    bvh_node(const hittable_list& list)
        : bvh_node(list.objects, 0, list.objects.size())
    {}

    // Constructor that takes a vector of hittable objects, a starting index and an ending index
    bvh_node(const std::vector<shared_ptr<hittable>>& src_objects, size_t start, size_t end);

    // Overridden hit method
    virtual bool hit(const Ray& r, double t_min, double t_max, hit_record& rec) const override;

    // Overridden bounding_box method
    virtual bool bounding_box(aabb& output_box) const override;

public: // left and right pointers to generic hittables (any primitive) allow us to split our hierarchy
    shared_ptr<hittable> left;
    shared_ptr<hittable> right;
    aabb box;
};

//implementation of bouding box method
bool bvh_node::bounding_box(aabb& output_box) const {
    output_box = box;
    return true;
}

//define hit method
bool bvh_node::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
    if (!box.hit(r, t_min, t_max))
        return false;

    bool hit_left = left->hit(r, t_min, t_max, rec);
    bool hit_right = right->hit(r, t_min, hit_left ? rec.t : t_max, rec);

    return hit_left || hit_right;
}

// Function to compare two hittable objects' bounding boxes
// based on the specified axis
inline bool box_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b, int axis) {
    aabb box_a;
    aabb box_b;

    if (!a->bounding_box(box_a) || !b->bounding_box(box_b))
        std::cerr << "No bounding box in bvh_node constructor.\n";

    return box_a.min()[axis] < box_b.min()[axis];
}

// Function to compare two hittable objects' bounding boxes along the x-axis
bool box_x_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return box_compare(a, b, 0);
}

// Function to compare two hittable objects' bounding boxes along the y-axis
bool box_y_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return box_compare(a, b, 1);
}

// Function to compare two hittable objects' bounding boxes along the z-axis
bool box_z_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return box_compare(a, b, 2);
}

// Constructor for BVH node which takes a vector of shared pointers to hittable objects, and start and end indices of the objects to be considered
bvh_node::bvh_node(const std::vector<shared_ptr<hittable>>& src_objects, size_t start, size_t end) {
    auto objects = src_objects; // Create a modifiable array of the source scene object

    // Choose a random axis to split on and select the appropriate comparator function
    int axis = random_int(0, 2);
    auto comparator = (axis == 0) ? box_x_compare
        : (axis == 1) ? box_y_compare
        : box_z_compare;

    size_t object_span = end - start;

    // If there is only one object, assign it to both left and right nodes
    if (object_span == 1) {
        left = right = objects[start];
    }
    // If there are two objects, sort them according to the selected comparator function and assign them to left and right nodes
    else if (object_span == 2) {
        if (comparator(objects[start], objects[start + 1])) {
            left = objects[start];
            right = objects[start + 1];
        }
        else {
            left = objects[start + 1];
            right = objects[start];
        }
    }
    // If there are more than two objects, sort them according to the selected comparator function and recursively build the left and right nodes
    else {
        std::sort(objects.begin() + start, objects.begin() + end, comparator);

        auto mid = start + object_span / 2;
        left = make_shared<bvh_node>(objects, start, mid);
        right = make_shared<bvh_node>(objects, mid, end);
    }

    // Compute the bounding boxes of the left and right nodes and compute the overall bounding box of the node
    aabb box_left, box_right;
    if (!left->bounding_box(box_left) || !right->bounding_box(box_right))
        std::cerr << "No bounding box in bvh_node constructor.\n";
    box = surrounding_box(box_left, box_right);
}