#pragma once
#include "geometry.h"
#include <chrono>

class vertexShader {
public:
	vertexShader() { //constructor to initialize the start time of the vertex shader
		t_start = std::chrono::high_resolution_clock::now();
	}
	void processVertex(Vec3f in, Vec3f* out) {
		auto passedTime = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t_start).count() / 1000; //calculate the time passed since the vertex shader started running
		//out->x = in.x +sinf(in.y*9 + passedTime)/9; //apply a sinusoidal function to the x-coordinate of the input vertex to make the model look groovy and dance haha
		out->x = in.x; //leave the x-coordinate of the input vertex unchanged
		out->y = in.y; //leave the y-coordinate of the input vertex unchanged
		out->z = in.z; //leave the z-coordinate of the input vertex unchanged
	}
	std::chrono::steady_clock::time_point t_start; //the start time of the vertex shader
};