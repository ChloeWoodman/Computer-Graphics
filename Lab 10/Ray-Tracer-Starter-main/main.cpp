#include "geometry.h"
#include "SDL.h"
#include "ray.h"
#include "hittable.h"
#include "common.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"
#include <fstream>
#include <chrono>

#define M_PI 3.14159265359

SDL_Window* window;
SDL_Renderer* renderer;
SDL_Surface* screen;

// Initializes the SDL library, creates a window, sets the renderer's draw color, and gets the screen surface.
void init() {
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window* window = SDL_CreateWindow(
        "Software Ray Tracer",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        640,
        480,
        0
    );

    screen = SDL_GetWindowSurface(window);

    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
}

// Sets the pixel at the specified (x, y) coordinate on the surface to the specified pixel value.
void putpixel(SDL_Surface* surface, int x, int y, Uint32 pixel)
{
    int bpp = surface->format->BytesPerPixel;
    /* Here p is the address to the pixel we want to set */
    Uint8* p = (Uint8*)surface->pixels + y * surface->pitch + x * bpp;

    switch (bpp) {
    case 1:
        *p = pixel;
        break;

    case 2:
        *(Uint16*)p = pixel;
        break;

    case 3:
        if (SDL_BYTEORDER == SDL_BIG_ENDIAN) {
            p[0] = (pixel >> 16) & 0xff;
            p[1] = (pixel >> 8) & 0xff;
            p[2] = pixel & 0xff;
        }
        else {
            p[0] = pixel & 0xff;
            p[1] = (pixel >> 8) & 0xff;
            p[2] = (pixel >> 16) & 0xff;
        }
        break;

    case 4:
        *(Uint32*)p = pixel;
        break;
    }
}

// method to ensure colours don’t go out of 8 bit range in RGB​
void clamp255(Vec3f& col) {
    col.x = (col.x > 255) ? 255 : (col.x < 0) ? 0 : col.x;
    col.y = (col.y > 255) ? 255 : (col.y < 0) ? 0 : col.y;
    col.z = (col.z > 255) ? 255 : (col.z < 0) ? 0 : col.z;
}

// Finds the t-value of the intersection of a ray with a sphere.
double hit_sphere(const Point3f& centre, double radius, const Ray& r) {
    Vec3f oc = r.origin() - centre;
    auto a = r.direction().dotProduct(r.direction());
    auto b = 2.0 * oc.dotProduct(r.direction());
    auto c = oc.dotProduct(oc) - radius * radius;
    auto discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return -1.0;
    }
    else {
        return (-b - sqrt(discriminant)) / (2.0 * a);
    }
}

// Determines the colour seen along a ray using recursive ray tracing.
Colour ray_colour(const Ray& r, const hittable& world, int depth) {
    hit_record rec;
    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0) return Colour(0, 0, 0);
    if (world.hit(r, 0.001, infinity, rec)) {  // bounded t from 0 --> inf
        Ray scattered;
        Colour attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_colour(scattered, world, depth - 1);
        return Colour(0, 0, 0);
    }

    Vec3f unit_direction = r.direction().normalize();
    auto t = 0.5 * (unit_direction.y + 1.0);
    return (1.0 - t) * Colour(1.0, 1.0, 1.0) + t * Colour(0.5, 0.7, 1.0) * 255;
}

// Initialize the program
int main(int argc, char** argv)
{
    init();
    // Image settings
    int spp = 10; // samples per pixel
    const auto aspect_ratio = 16.0 / 9.0; // width to height ratio of the image
    const int image_width = screen->w; // width of the output image
    const int image_height = static_cast<int>(image_width / aspect_ratio); // height of the output image
    const float scale = 1.0f / spp; // scale the colour values according to the number of samples per pixel
    const int max_depth = 50; // maximum depth for recursive ray tracing

    // Camera settings
    auto viewport_height = 2.0; // height of the viewport
    auto viewport_width = aspect_ratio * viewport_height; // width of the viewport
    auto focal_length = 1.0; // focal length of the camera
    auto origin = Point3f(0, 0, 0); // position of the camera
    auto horizontal = Vec3f(viewport_width, 0, 0); // horizontal vector of the camera viewport
    auto vertical = Vec3f(0, viewport_height, 0); // vertical vector of the camera viewport
    auto lower_left_corner = origin - horizontal / 2 - vertical / 2 - Vec3f(0, 0, focal_length); // lower left corner of the camera viewport
    camera cam;

    // Colour values for the scene
    const Colour white(255, 255, 255);
    const Colour black(0, 0, 0);
    const Colour red(255, 0, 0);

    // Create the scene
    hittable_list world;
    auto material_ground = make_shared<lambertian>(Colour(0.8, 0.8, 0.0));
    auto material_center = make_shared<lambertian>(Colour(0.7, 0.3, 0.3));
    auto material_left = make_shared<metal>(Colour(0.8, 0.8, 0.8), 0.3);
    auto material_right = make_shared<metal>(Colour(0.8, 0.6, 0.2), 1.0);

    // Add objects to the world
    world.add(make_shared<Sphere>(Point3f(0.0, -100.5, -1.0), 100.0, material_ground));
    world.add(make_shared<Sphere>(Point3f(0.0, 0.0, -1.0), 0.5, material_center));
    world.add(make_shared<Sphere>(Point3f(-1.0, 0.0, -1.0), 0.5, material_left));
    world.add(make_shared<Sphere>(Point3f(1.0, 0.0, -1.0), 0.5, material_right));

    double t;
    Colour pix_col(black);

    // Render loop
    SDL_Event e;
    bool running = true;
    while (running) {

        auto t_start = std::chrono::high_resolution_clock::now();

        // Clear the back buffer, pixel data on surface and depth buffer
        SDL_FillRect(screen, nullptr, SDL_MapRGB(screen->format, 0, 0, 0));
        SDL_RenderClear(renderer);

        // The following loop iterates over every pixel in the screen and renders the corresponding ray traced colour.
        // Loop through each pixel in the image
        for (int y = screen->h - 1; y >= 0; --y) { // start from the top left

            //std::cerr << "\rScanlines remaining: " << y << std::flush;
            for (int x = 0; x < screen->w; ++x) {
                // Reset the colour for every pixel.
                pix_col = black;

                // Perform spp samples at each pixel to accumulate colours.
                for (int s = 0; s < spp; s++) {
                    auto u = double(x + random_double()) / (image_width - 1);
                    auto v = double(y + random_double()) / (image_height - 1);
                    Ray ray = cam.get_ray(u, v);

                    // Accumulate colours for every sample.
                    pix_col = pix_col + ray_colour(ray, world, max_depth);
                }

                // Scale cumulative colour values accordingly to spp.
                pix_col /= 255.f * spp;

                // Apply gamma correction.
                pix_col.x = sqrt(pix_col.x);
                pix_col.y = sqrt(pix_col.y);
                pix_col.z = sqrt(pix_col.z);

                // Rescale pixel colour values to the range 0-255.
                pix_col *= 255;

                // Map the colour to an SDL-compatible format.
                Uint32 colour = SDL_MapRGB(screen->format, pix_col.x, pix_col.y, pix_col.z);

                // Put the pixel on the screen.
                putpixel(screen, x, y, colour);
            }
        }

        // Calculate the frame render time and print it out
        auto t_end = std::chrono::high_resolution_clock::now();
        auto passedTime = std::chrono::duration<double, std::milli>(t_end - t_start).count();
        std::cerr << "Frame render time:  " << passedTime << " ms" << std::endl;

        // Create an SDL texture from the rendered image.
        SDL_Texture* texture = SDL_CreateTextureFromSurface(renderer, screen);
        if (texture == NULL) {
            fprintf(stderr, "CreateTextureFromSurface failed: %s\n", SDL_GetError());
            exit(1);
        }

        // Free the memory associated with the screen surface.
        SDL_FreeSurface(screen);

        // Render the texture to the screen.
        SDL_RenderCopyEx(renderer, texture, NULL, NULL, 0, 0, SDL_FLIP_VERTICAL);
        SDL_RenderPresent(renderer);

        // Destroy the texture.
        SDL_DestroyTexture(texture);

        if (SDL_PollEvent(&e))
        {
            switch (e.type) {
            case SDL_KEYDOWN:
                switch (e.key.keysym.sym) {
                case SDLK_ESCAPE:
                    running = false;
                    break;
                }
                break;
            }
        }
    }
    return 0;
}