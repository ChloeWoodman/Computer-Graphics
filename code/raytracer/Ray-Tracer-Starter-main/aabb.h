#pragma once
#include "common.h"

// Definition of aabb class
class aabb {
public:
    // Default constructor
    aabb() {}
    aabb(const Point3f& mini, const Point3f& maxi) { minimum = mini; maximum = maxi; } // Constructor that sets the minimum and maximum points.

    Point3f min() const { return minimum; } // Returns the minimum point of the AABB.
    Point3f max() const { return maximum; } // Returns the maximum point of the AABB.

    bool hit(const Ray& r, double t_min, double t_max) const { // Determines if a ray hits the AABB.
        for (int a = 0; a < 3; a++) { // For each axis:
            auto t0 = fmin((minimum[a] - r.origin()[a]) / r.direction()[a], // Calculate the intersection with the near face.
                (maximum[a] - r.origin()[a]) / r.direction()[a]); // Calculate the intersection with the far face.
            auto t1 = fmax((minimum[a] - r.origin()[a]) / r.direction()[a],
                (maximum[a] - r.origin()[a]) / r.direction()[a]);
            t_min = fmax(t0, t_min); // Update t_min to the maximum of the near face intersection and the previous value of t_min.
            t_max = fmin(t1, t_max); // Update t_max to the minimum of the far face intersection and the previous value of t_max.
            if (t_max <= t_min) // If there is no overlap between the intersection intervals, the ray misses the AABB.
                return false;
        }
        return true; // If the ray intersects the AABB on all three axes, return true.
    }

    Point3f minimum; // The minimum point of the AABB.
    Point3f maximum; // The maximum point of the AABB.
};

aabb surrounding_box(aabb box0, aabb box1) { // Computes the surrounding AABB of two given boxes.
    Point3f small(fmin(box0.min().x - 1e-3, box1.min().x - 1e-3), // Computes the minimum point of the surrounding AABB.
        fmin(box0.min().y - 1e-3, box1.min().y - 1e-3),
        fmin(box0.min().z - 1e-3, box1.min().z - 1e-3));

    Point3f big(fmax(box0.max().x + 1e-3, box1.max().x + 1e-3), // Computes the maximum point of the surrounding AABB.
        fmax(box0.max().y + 1e-3, box1.max().y + 1e-3),
        fmax(box0.max().z + 1e-3, box1.max().z + 1e-3));

    return aabb(small, big); // Constructs and returns the surrounding AABB.
}
