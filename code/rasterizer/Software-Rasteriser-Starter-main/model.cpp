#include "model.h"

Model::Model(const char* filename) : verts_(), faces_(), vts_(), vns_() {
    std::ifstream in;
    in.open(filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;

        // Read in vertex data
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i = 0; i < 3; i++) iss >> v[i];
            verts_.push_back(v);
        }
        // Read in texture coordinate data
        else if (!line.compare(0, 3, "vt ")) {
            iss >> trash;
            iss >> trash;
            Vec2f vt;
            for (int i = 0; i < 2; i++) iss >> vt[i];
            vts_.push_back(vt);
        }

        // Read in vertex normal data
        else if (!line.compare(0, 3, "vn ")) {
            iss >> trash;
            iss >> trash;
            Vec3f vn;
            for (int i = 0; i < 3; i++) iss >> vn[i];
            vns_.push_back(vn); //holds values for vertex normals
        }

        // Read in face data
        else if (!line.compare(0, 2, "f ")) {
            std::vector<int> f;
            std::vector<int> vns;
            std::vector<int> uvs;
            int itrash, v_idx, vt_idx, vn_idx;
            iss >> trash;
            while (iss >> v_idx >> trash >> vt_idx >> trash >> vn_idx) { 
                v_idx--, vt_idx--, vn_idx--; //in wavefront obj the indices start at 1 instead of 0
                f.push_back(v_idx);
                vns.push_back(vn_idx);
                uvs.push_back(vt_idx);
            }
            faces_.push_back(f);
            vnorms_.push_back(vns);
            uvs_.push_back(uvs);
        }    
    }

    // Print out information about the loaded model
    std::cerr << "- v# " << verts_.size() << " f# " << faces_.size() << " vn# " << vns_.size() << std::endl;
}

Model::~Model() {
}

int Model::nverts() {
    return (int)verts_.size();
}

int Model::nnorms() {
    return (int)vns_.size();
}

int Model::nfaces() {
    return (int)faces_.size();
}

std::vector<int> Model::face(int idx) {
    return faces_[idx];
}

std::vector<int> Model::vertex_normal(int idx) {
    return vnorms_[idx];
}

std::vector<int> Model::uv_coord(int idx) {
    return uvs_[idx];
}

std::vector<int> Model::vt_indices(int idx) {
    return uvs_[idx];
}



Vec3f Model::vert(int i) {
    return verts_[i];
}

Vec2f Model::vt(int i) {
    return vts_[i];
}

Vec3f Model::vn(int i) {
    return vns_[i];
}