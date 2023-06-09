#pragma once
#include <cmath> // Required for mathematical functions
#include <limits> // Required for numeric_limits
#include <memory> // Required for smart pointers
#include <cstdlib> // Required for rand function

//using declarations for convenience and readability
using std::shared_ptr;
using std::make_shared;
using std::sqrt;

//constants used in the program
const double infinity = std::numeric_limits<double>::infinity(); // A constant representing positive infinity
const double pi = 3.1415926535897932385; // A constant representing the value of pi

//utility functions used in the program
inline double degrees_to_radians(double degrees) { // Converts degrees to radians
	return degrees * pi / 180.0;
}
inline double random_double() { // Returns a random real in [0,1]
	return rand() / (RAND_MAX + 1.0);
}
inline double random_double(double min, double max) { // Returns a random real in [min, max]
	return min + (max - min) * random_double();
}

#include <random>

float random_float(float min, float max) {
	static std::uniform_real_distribution<float> distribution(min, max);
	static std::mt19937 generator;
	return distribution(generator);
}


// Helper methods for the bvh_node constructor
//creating a generic comparatot that will return true if argument 1 < 1 second
//..given an additional axis index argument
//This then defined the axis-specific comparison functions that utilize the generic comparison function
inline int random_int(int min, int max) {
	// Returns a random integer in [min,max].
	return static_cast<int>(random_double(min, max + 1));
}


// This function clamps a value between a minimum and maximum value
// If the value is less than the minimum, it returns the minimum
// If the value is greater than the maximum, it returns the maximum
// Otherwise, it returns the value itself
inline double clamp(double x, double min, double max) {
	if (x < min) return min;
	if (x > max) return max;
	return x;
}


//common header files used in the program
#include "ray.h" // Contains the declaration of the Ray class
#include "geometry.h" // Contains the declaration of some vector and point classes