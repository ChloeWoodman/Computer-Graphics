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

//common header files used in the program
#include "ray.h" // Contains the declaration of the Ray class
#include "geometry.h" // Contains the declaration of some vector and point classes