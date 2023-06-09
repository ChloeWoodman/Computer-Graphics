#include <vector>
#include "geometry.h"

class Model {
private:
    std::vector<Vec3f> verts_;
    std::vector<std::vector<int>> faces_;
    std::vector<Vec2f> vts_;
    std::vector<std::vector<int>> vt_indices_; // Vector of vectors to store texture coordinate indices for each face

public:
    Model(const char* filename);
    ~Model();
    int nverts();
    int nfaces();
    Vec3f vert(int i);
    Vec2f vt(int i);
    std::vector<int> face(int idx);
    std::vector<int> vt_indices(int idx); // Function to return the texture coordinate indices for the face at index idx
};
