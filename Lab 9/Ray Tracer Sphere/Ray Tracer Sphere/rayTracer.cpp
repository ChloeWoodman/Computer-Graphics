#include "rayTracer.h"
#include <fstream>
#include <cmath>

// Define a 3D vector class
struct Vec3 {
	double x, y, z;
	Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
	// Define vector operations
	Vec3 operator + (const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
	Vec3 operator - (const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
	Vec3 operator * (double d) const { return Vec3(x * d, y * d, z * d); }
	Vec3 operator / (double d) const { return Vec3(x / d, y / d, z / d); }
	// Normalize the vector
	Vec3 normalize() const {
		double mg = sqrt(x * x + y * y + z * z);
		return Vec3(x / mg, y / mg, z / mg);
	}
};

// Calculate the dot product of two vectors
double dot(const Vec3& a, const Vec3& b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Define a ray class
struct Ray {
	Vec3 o, d;
	Ray(const Vec3& o, const Vec3& d) : o(o), d(d) {}
};

// Define a sphere class
struct Sphere {
	Vec3 c;
	double r;
	Sphere(const Vec3& c, double r) : c(c), r(r) {}
	// Calculate the normal vector at a point on the sphere
	Vec3 getNormal(const Vec3& pi) const { return (pi - c) / r; }
	// Check if a ray intersects the sphere and calculate the intersection point
	bool intersect(const Ray& ray, double& t) const {
		const Vec3 o = ray.o;
		const Vec3 d = ray.d;
		const Vec3 oc = o - c;
		const double b = 2 * dot(oc, d);
		const double a = dot(d, d);
		const double c = dot(oc, oc) - r * r;
		double disc = b * b - 4 * a * c;
		if (disc < 1e-4) return false;
		disc = sqrt(disc);
		const double t0 = (-b - disc) / (2 * a);
		const double t1 = (-b + disc) / (2 * a);
		t = (t0 < 1) ? t0 : t1;
		return true;
	}
};

// Clamp the color values to the range [0, 255]
void clamp255(Vec3& col) {
	col.x = (col.x > 255) ? 255 : (col.x < 0) ? 0 : col.x;
	col.y = (col.y > 255) ? 255 : (col.y < 0) ? 0 : col.y;
	col.z = (col.z > 255) ? 255 : (col.z < 0) ? 0 : col.z;
}

int main() {
	const int H = 500;
	const int W = 500;

	const Vec3 white(255, 255, 255);
	const Vec3 black(0, 0, 0);
	const Vec3 red(255, 0, 0);

	const Sphere sphere(Vec3(W * 0.5, H * 0.5, 50), 50);
	const Sphere light(Vec3(0, 0, 50), 1);

	std::ofstream out("out.ppm");
	out << "P3\n" << W << ' ' << H << ' ' << "255\n";

	double t;
	Vec3 pix_col(black);

	// Loop through each pixel in the image
	for (int y = 0; y < H; ++y) {
		for (int x = 0; x < W; ++x) {
			pix_col = black;

			const Ray ray(Vec3(x, y, 0), Vec3(0, 0, 1));
			// Check if the ray intersects the sphere
			if (sphere.intersect(ray, t)) {
				const Vec3 pi = ray.o + ray.d * t;
				const Vec3 L = light.c - pi;
				const Vec3 N = sphere.getNormal(pi);
				const double dt = dot(L.normalize(), N.normalize());

				// Calculate the pixel color based on the intersection
				pix_col = (red + white * dt) * 0.5;
				clamp255(pix_col);
			}
			// Write the pixel color to the output file
			out << (int)pix_col.x << ' '
				<< (int)pix_col.y << ' '
				<< (int)pix_col.z << '\n';
		}
	}
}