#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
Model *model = NULL;


const int width  = 800;
const int height = 800;


// Write the line method here

void plotLine(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor colour) {
    float y;
    float t;
    bool steep = false;
    //check if the line is steep and swap the x and y values (transpose calculations)
    if ((abs(x0-x1)) < (abs(y0-y1))) {
        std::swap(x0, y0);
        std::swap(x1, y1);

        steep = true;
    }

    //make the line left to right
    if (x0 > x1){
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    for (int x = x0; x <= x1; x++) {
        t = (x-(float)x0) / ((float)x1 - (float)x0);
        y = y0 * (1 - t) + (y1*t);
        //if the image is transposed, de-transpose it
        if (steep){
            image.set(y, x, colour);
        }
        else {
            image.set(x, y, colour);
        }
        
    }
}

int main(int argc, char** argv) {

    TGAImage image(width, height, TGAImage::RGB);

    if (2 == argc) {
        model = new Model(argv[1]);
    }
    else {
        model = new Model("cc.obj");
    }

    for (int i = 0; i < model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        for (int j = 0; j < 3; j++) {
            Vec3f v0 = model->vert(face[j]);
            Vec3f v1 = model->vert(face[(j + 1) % 3]);
            int x0 = (v0.x + 1.) * width / 2.;
            int y0 = (v0.y + 1.) * height / 2.;
            int x1 = (v1.x + 1.) * width / 2.;
            int y1 = (v1.y + 1.) * height / 2.;
            plotLine(x0, y0, x1, y1, image, white);
        }
    }

    image.flip_vertically(); // we want to have the origin at the left bottom corner of the image
    //plotLine(13, 20, 80, 40, image, white);
    //plotLine(20, 13, 40, 80, image, red);
    //plotLine(80, 40, 13, 20, image, red);
    image.write_tga_file("output.tga");
    delete model;
    return 0;
}
