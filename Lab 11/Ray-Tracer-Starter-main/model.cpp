#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

Model::Model(const char* filename) : verts_(), faces_(), vts_(), vt_indices_() {
    std::ifstream in;
    in.open(filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i = 0; i < 3; i++) iss >> v[i];
            verts_.push_back(v);
        }
        else if (!line.compare(0, 3, "vt ")) {
            iss >> trash;
            Vec2f vt;
            for (int i = 0; i < 2; i++) iss >> vt[i];
            vts_.push_back(vt);
        }
        else if (!line.compare(0, 2, "f ")) {
            std::vector<int> f;
            std::vector<int> vt_f; // Vector to store texture coordinate indices for the current face
            int itrash, idx, vt_idx; // New variable vt_idx to store texture coordinate index
            iss >> trash;
            while (iss >> idx >> trash >> vt_idx >> trash >> itrash) { // Read in v_i to idx, vt_i to vt_idx, and discard /vn_i
                idx--;
                vt_idx--; // Texture coordinate indices also start at 1, so we need to subtract 1 to start at 0
                f.push_back(idx);
                vt_f.push_back(vt_idx); // Add the texture coordinate index to the vt_f vector
            }
            faces_.push_back(f);
            vt_indices_.push_back(vt_f); // Add the vt_f vector to the vt_indices_ vector
        }
    }
    std::cerr << "# v# " << verts_.size() << " f# " << faces_.size() << " vt# " << vts_.size() << std::endl;
}

Model::~Model() {
}

int Model::nverts() {
    return (int)verts_.size();
}

int Model::nfaces() {
    return (int)faces_.size();
}

std::vector<int> Model::face(int idx) {
    return faces_[idx];
}

Vec3f Model::vert(int i) {
    return verts_[i];
}

Vec2f Model::vt(int i) {
    return vts_[i];
}