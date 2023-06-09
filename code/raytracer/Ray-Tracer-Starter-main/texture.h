#pragma once
// stb from https://github.com/nothings/stb.git
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "perlin.h"
#include <iostream>
#include <cmath>

// This class defines a base interface for all textures.
class texture {
public:
    // This function returns the color of the texture at the given (u, v) texture coordinates and the point p on the surface.
    virtual Colour value(double u, double v, const Point3f& p) const = 0;
};

class noise_texture : public texture {
public:
    noise_texture() {}
    noise_texture(double sc) : scale(sc) {}

    virtual Colour value(double u, double v, const Point3f& p) const override {
        return Colour(1, 1, 1) * 0.5 * (1 + sin(scale * p.z + 10 * noise.turb(p)));
    }

public:
    perlin noise;
    double scale;
};

// This class defines a texture that has a solid color.
class solid_color : public texture {
public:
    solid_color() {}
    solid_color(Colour c) : colour_value(c) {}
    // This constructor creates a solid color texture from three color values (red, green, and blue).
    solid_color(double red, double green, double blue)
        : solid_color(Colour(red, green, blue)) {}

    // This function returns the color of the solid color texture at the given (u, v) texture coordinates and the point p on the surface.
    virtual Colour value(double u, double v, const Vec3f& p) const override {
        return colour_value;
    }

private:
    Colour colour_value;
};

// This class defines a checkerboard texture.
class checker_texture : public texture {
public:
    checker_texture() {}

    // This constructor creates a checkerboard texture from two other textures (even and odd).
    checker_texture(shared_ptr<texture> _even, shared_ptr<texture> _odd)
        : even(_even), odd(_odd) {}

    // This constructor creates a checkerboard texture from two colors (c1 and c2).
    checker_texture(Colour c1, Colour c2)
        : even(make_shared<solid_color>(c1)), odd(make_shared<solid_color>(c2)) {}

    // This function returns the color of the checkerboard texture at the given (u, v) texture coordinates and the point p on the surface.
    virtual Colour value(double u, double v, const Point3f& p) const override {
        auto sines = sin(10 * p.x) * sin(10 * p.y) * sin(10 * p.z);
        if (sines < 0)
            return odd->value(u, v, p);
        else
            return even->value(u, v, p);
    }

public:
    shared_ptr<texture> odd;
    shared_ptr<texture> even;
};

// This class defines an image texture.
class image_texture : public texture {
public:
    const static int bytes_per_pixel = 3;

    image_texture()
        : data(nullptr), width(0), height(0), bytes_per_scanline(0) {}

    // This constructor creates an image texture from a file.
    image_texture(const char* filename) {
        auto components_per_pixel = bytes_per_pixel;

        data = stbi_load(
            filename, &width, &height, &components_per_pixel, components_per_pixel);

        if (!data) {
            std::cerr << "ERROR: Could not load texture image file '" << filename << "'.\n";
            width = height = 0;
        }

        bytes_per_scanline = bytes_per_pixel * width;
    }

    ~image_texture() {
        delete data;
    }

    // This function returns the color of the image texture at the given (u, v) texture coordinates and the point p on the surface.
    virtual Colour value(double u, double v, const Vec3f& p) const override {
        // If we have no texture data, then return solid cyan as a debugging aid.
        if (data == nullptr)
            return Colour(0, 1, 1);

        // Clamp input texture coordinates to [0,1] x [1,0]
        u = clamp(u, 0.0, 1.0);
        v = 1.0 - clamp(v, 0.0, 1.0);  // Flip V to image coordinates

        auto i = static_cast<int>(u * width);
        auto j = static_cast<int>(v * height);

        // Clamp integer mapping, since actual coordinates should be less than 1.0
        if (i >= width)  i = width - 1;
        if (j >= height) j = height - 1;

        const auto color_scale = 1.0 / 255.0;
        auto pixel = data + j * bytes_per_scanline + i * bytes_per_pixel;

        return Colour(color_scale * pixel[0], color_scale * pixel[1], color_scale * pixel[2]);
    }

private:
    unsigned char* data;
    int width, height;
    int bytes_per_scanline;
};