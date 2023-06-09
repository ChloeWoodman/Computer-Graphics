#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "geometry.h"

class Model {
private:
    std::vector<Vec3f> verts_;
    std::vector<std::vector<int>> faces_;
    std::vector<Vec2f> vts_;
    std::vector<std::vector<int>> uvs_;
    std::vector<Vec3f> vns_;
    std::vector<std::vector<int>> vnorms_;

public:
    Model(const char* filename);
    ~Model();
    int nverts();
    int nnorms();
    int nfaces();
    Vec3f vert(int i);
    Vec2f vt(int i);
    Vec3f vn(int i);
    std::vector<int> face(int idx);
    std::vector<int> vertex_normal(int idx);
    std::vector<int> uv_coord(int idx);
    std::vector<int> vt_indices(int idx);
};