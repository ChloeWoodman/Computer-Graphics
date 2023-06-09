#include "common.h"
#include <cmath>

// Class definition for the perlin noise generator
class perlin {
public:
    // Constructor
    perlin() {
        // Initialization of random vectors at the corners of the unit cube
        ranvec = new Vec3f[point_count];
        for (int i = 0; i < point_count; ++i) {
            ranvec[i] = (Vec3f::random(-1, 1)).normalize();
        }

        // Initialization of permutation arrays to randomize access to the vectors
        perm_x = perlin_generate_perm();
        perm_y = perlin_generate_perm();
        perm_z = perlin_generate_perm();
    }

    // Destructor
    ~perlin() {
        // Deallocating memory
        delete[] ranvec;
        delete[] perm_x;
        delete[] perm_y;
        delete[] perm_z;
    }

    // Perlin noise function
    double noise(const Point3f& p) const {
        // Computing the unit cube that contains the input point
        auto u = p.x - floor(p.x);
        auto v = p.y - floor(p.y);
        auto w = p.z - floor(p.z);
        auto i = static_cast<int>(floor(p.x));
        auto j = static_cast<int>(floor(p.y));
        auto k = static_cast<int>(floor(p.z));
        Vec3f c[2][2][2];

        // Fetching the eight gradient vectors for the corners of the cube
        for (int di = 0; di < 2; di++) {
            for (int dj = 0; dj < 2; dj++) {
                for (int dk = 0; dk < 2; dk++) {
                    c[di][dj][dk] = ranvec[
                        perm_x[(i + di) & 255] ^
                            perm_y[(j + dj) & 255] ^
                            perm_z[(k + dk) & 255]
                    ];
                }
            }
        }

        // Interpolating the eight gradient values to get the final noise value
        return perlin_interp(c, u, v, w);
    }

    // Turbulence function, which is a fractal sum of noise
    double turb(const Point3f& p, int depth = 7) const {
        auto accum = 0.0;
        auto temp_p = p;
        auto weight = 1.0;

        for (int i = 0; i < depth; i++) {
            accum += weight * noise(temp_p);
            weight *= 0.5; // Decreasing the weight by half each iteration
            temp_p *= 2; // Increasing the frequency by two each iteration
        }
        return fabs(accum); // Taking the absolute value to make the function symmetric about the origin
    }

private:
    static const int point_count = 256;
    Vec3f* ranvec;
    int* perm_x;
    int* perm_y;
    int* perm_z;

    // Method for generating a random permutation of {0,1,...,point_count-1}
    static int* perlin_generate_perm() {
        auto p = new int[point_count];

        for (int i = 0; i < perlin::point_count; i++)
            p[i] = i;

        permute(p, point_count);

        return p;
    }

    // Method for permuting an array
    static void permute(int* p, int n) {
        for (int i = n - 1; i > 0; i--) {
            int target = random_int(0, i);
            int tmp = p[i];
            p[i] = p[target];
            p[target] = tmp;
        }
    }

    // Method for performing trilinear interpolation of gradient values
    static double perlin_interp(Vec3f c[2][2][2], double u, double v, double w) {
        auto uu = u * u * (3 - 2 * u);
        auto vv = v * v * (3 - 2 * v);
        auto ww = w * w * (3 - 2 * w);
        auto accum = 0.0;

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 2; k++) {
                    Vec3f weight_v(u - i, v - j, w - k);
                    // Trilinear interpolation formula
                    accum += (i * uu + (1 - i) * (1 - uu))
                        * (j * vv + (1 - j) * (1 - vv))
                        * (k * ww + (1 - k) * (1 - ww))
                        * c[i][j][k].dotProduct(weight_v);
                }
            }
        }
        return accum;
    }
};
