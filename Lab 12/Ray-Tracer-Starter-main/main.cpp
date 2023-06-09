#include "geometry.h"
#include "SDL.h"
#include "ray.h"
#include "hittable.h"
#include "common.h"
#include "hittable_list.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"
#include "threadpool.h"
#include "triangle.h"
#include "model.h"
#include "bvh.h"
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
Colour ray_colour(const Ray& r, const Colour& background, const hittable& world, int depth) {
    hit_record rec;
    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0) return Colour(0, 0, 0);
    if (!world.hit(r, 0.001, infinity, rec)) return background; // bounded t from 0 --> inf

    Ray scattered;
    Colour attenuation;
    Colour emitted = rec.mat_ptr->emitted();

    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
        return emitted;

    return emitted + attenuation * ray_colour(scattered, background, world, depth - 1);
}


    /*Vec3f unit_direction = r.direction().normalize();
    auto t = 0.5 * (unit_direction.y + 1.0);
    return (1.0 - t) * Colour(1.0, 1.0, 1.0) + t * Colour(0.5, 0.7, 1.0) * 255;*/
//}

// Test rendering a new scene
hittable_list random_scene() {
    hittable_list world;
    auto ground_material = make_shared<lambertian>(Colour(0.5, 0.5, 0.5));
    world.add(make_shared<Sphere>(Point3f(0, -1000, 0), 1000, ground_material));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            Point3f centre(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());
            if ((centre - Point3f(4, 0.2, 0)).length() > 0.9) {
                shared_ptr<material> sphere_material;
                if (choose_mat < 0.8) { // diffuse
                    auto albedo = Colour::random() * Colour::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<Sphere>(centre, 0.2, sphere_material));
                }
                else if (choose_mat < 0.95) { // metal
                    auto albedo = Colour::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<Sphere>(centre, 0.2, sphere_material));
                }
                else { // glass  
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<Sphere>(centre, 0.2, sphere_material));
                }
            }
        }
    }
    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<Sphere>(Point3f(0, 1, 0), 1.0, material1));
    auto material2 = make_shared<lambertian>(Colour(0.4, 0.2, 0.1));
    world.add(make_shared<Sphere>(Point3f(-4, 1, 0), 1.0, material2));
    auto material3 = make_shared<metal>(Colour(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<Sphere>(Point3f(4, 1, 0), 1.0, material3));
    auto material4 = make_shared<diffuse_light>(Colour(255, 255, 255));
    world.add(make_shared<Sphere>(Point3f(0, 5, 0), 1.0, material4));
    
    //return world;                                       // without bvh
    return hittable_list(make_shared<bvh_node>(world)); // with bvh
}

void LineRender(SDL_Surface* screen, hittable_list world, int y, int spp, int max_depth, camera* cam) {
    const float aspect_ratio = 16.0 / 9;
    const int image_width = screen->w;
    const int image_height = static_cast<int>(image_width / aspect_ratio);

    const Colour black(0, 0, 0);
    Colour pix_col(black);

    Colour background(0, 0, 0);
    for (int x = 0; x < screen->w; ++x) {
        pix_col = black;
        for (int s = 0; s < spp; s++) {
            auto u = double(x + random_double()) / (image_width - 1);
            auto v = double(y + random_double()) / (image_height - 1);
            Ray ray = cam->get_ray(u, v);
            Vec3f unit_direction = ray.direction().normalize();
            auto t = 0.5 * (unit_direction.y + 1.0);
            background = (1.0 - t) * Colour(1.0, 1.0, 1.0) + t * Colour(0.5, 0.7, 1.0) * 255;
            pix_col = pix_col + ray_colour(ray, background, world, max_depth);
        }
        pix_col /= 255.f * spp;
        pix_col.x = sqrt(pix_col.x);
        pix_col.y = sqrt(pix_col.y);
        pix_col.z = sqrt(pix_col.z);
        pix_col *= 255;
        Uint32 colour = SDL_MapRGB(screen->format, pix_col.x, pix_col.y, pix_col.z);
        putpixel(screen, x, y, colour);
    }
}

hittable_list test_scene() {
    hittable_list world;

    //Use your own or the model class from the rasterisation lab
    Model* model = new Model("cc_t.obj");

    Vec3f transform(0, 0.8, 0);
    auto glass = make_shared<dielectric>(1.5);
    for (uint32_t i = 0; i < model->nfaces(); ++i) {
        const Vec3f& v0 = model->vert(model->face(i)[0]);
        const Vec3f& v1 = model->vert(model->face(i)[1]);
        const Vec3f& v2 = model->vert(model->face(i)[2]);
        world.add(make_shared<triangle>(v0 + transform, v1 + transform, v2 + transform, glass));
    }

    transform = Vec3f(1.2, 0.8, 0);
    auto mat_diffuse = make_shared<lambertian>(Colour(0.4, 0.2, 0.1));
    for (uint32_t i = 0; i < model->nfaces(); ++i) {
        const Vec3f& v0 = model->vert(model->face(i)[0]);
        const Vec3f& v1 = model->vert(model->face(i)[1]);
        const Vec3f& v2 = model->vert(model->face(i)[2]);
        world.add(make_shared<triangle>(v0 + transform, v1 + transform, v2 + transform, mat_diffuse));
    }

    transform = Vec3f(-1.2, 0.8, 0);
    auto mat_metal = make_shared<metal>(Colour(0.7, 0.6, 0.5), 0.0);
    mat_diffuse = make_shared<lambertian>(Colour(0.8, 0.4, 0.7));
    for (uint32_t i = 0; i < model->nfaces(); ++i) {
        const Vec3f& v0 = model->vert(model->face(i)[0]);
        const Vec3f& v1 = model->vert(model->face(i)[1]);
        const Vec3f& v2 = model->vert(model->face(i)[2]);
        world.add(make_shared<triangle>(v0 + transform, v1 + transform, v2 + transform, mat_metal));
    }

    auto ground_material = make_shared<lambertian>(Colour(0.5, 0.5, 0.5));
    world.add(make_shared<Sphere>(Point3f(0, -1000, 0), 1000, ground_material));
    auto material4 = make_shared<diffuse_light>(Colour(255, 255, 255));
    world.add(make_shared<Sphere>(Point3f(0, 5, 0), 1.0, material4));
   
    //return world;                                       // without bvh
    return hittable_list(make_shared<bvh_node>(world)); // with bvh
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

    // Colour values for the scene
    const Colour white(255, 255, 255);
    const Colour black(0, 0, 0);
    const Colour red(255, 0, 0);

    // Create the scene
    //auto world = random_scene();
    auto world = test_scene();

    // Camera settings
    Point3f lookfrom(13, 5, 23);
    Point3f lookat(0, 0, 0);
    Vec3f vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.15;
    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);
    
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

        //multithreading
        {
            t_start = std::chrono::high_resolution_clock::now();
            ThreadPool pool(std::thread::hardware_concurrency());

            int start = screen->h - 1;
            int step = screen->h / std::thread::hardware_concurrency();
            for (int y = 0; y < screen->h - 1; y++)
            {
                //pool.Enqueue(std::bind(screen, world, y, spp, max_depth, &cam));
                pool.Enqueue(std::bind(LineRender, screen, world, y, spp, max_depth, &cam));
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