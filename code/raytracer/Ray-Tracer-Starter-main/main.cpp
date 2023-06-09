#include "geometry.h"
#include "SDL.h"
#include "SDL_mixer.h"
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
#include "moving_sphere.h"
#include "moving_triangle.h"
#include "aarect.h"
#include "box.h"
#include "constant_medium.h"
#include <fstream>
#include <chrono>

#define M_PI 3.14159265359

SDL_Window* window;
SDL_Renderer* renderer;
SDL_Surface* screen;

// Colour values for the scene
const Colour white(255, 255, 255);
const Colour black(0, 0, 0);
const Colour red(255, 0, 0);

hittable_list world;

// Camera settings
auto aperture = 0.f;
Vec3f vup(0, 1, 0);
auto vfov = 20.0;
auto dist_to_focus = 10.0;
Colour background(0, 0, 0);
float camAngleX = 0.f, camAngleY = 0.f;
float camX = 0.0f, camY = 0.0f, camZ = 0.0f;
float distance = 3.0f; //straight line distance between camera and look at point
Point3f lookfrom = Point3f(camX + -180.122f, camY + 133.769f, camZ + 202.863f);
Point3f lookat = Point3f(-140, 150, 0);

// Initializes the SDL library, creates a window, sets the renderer's draw color, and gets the screen surface.
void init() {
	SDL_Init(SDL_INIT_VIDEO);

	SDL_Window* window = SDL_CreateWindow(
		"Software Ray Tracer",
		SDL_WINDOWPOS_UNDEFINED,
		SDL_WINDOWPOS_UNDEFINED,
		1920,
		1080,
		0
	);

	screen = SDL_GetWindowSurface(window);

	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
}

// Loads the music file and plays it
void playMusic() {
	// Initialize SDL_mixer
	if (Mix_OpenAudio(44100, MIX_DEFAULT_FORMAT, 2, 2048) < 0) {
		std::cout << "could not intialise mixer!" << std::endl;
		std::cout << SDL_GetError() << std::endl;
		return;
	}

	// Load music file
	//using audio from this music file https://www.youtube.com/watch?v=-5ySelLIAgQ
	Mix_Music* music = Mix_LoadMUS("./audio/Alone_With_Your_Voices.wav");
	if (!music) {
		// Handle music loading error
		printf("Failed to load music: %s\n", Mix_GetError());
		return;
	}

	// Play music
	if (Mix_PlayMusic(music, -1) == -1) {
		// Handle music playing error
		printf("Failed to play music: %s\n", Mix_GetError());
		return;
	}

	// Music is playing
	printf("Music is playing...\n");
}


// Sets the pixel at the specified (x, y) coordinate on the surface to the specified pixel value.
void putpixel(SDL_Surface* surface, int x, int y, Uint32 pixel)
{
	int bpp = surface->format->BytesPerPixel;
	/* Here p is the address to the pixel we want to set*/
	Uint8* p = (Uint8*)surface->pixels + y * surface->pitch + x * bpp;
	switch (bpp) {
	case 1:
		*p = pixel; // Set pixel for 8-bit surface
		break;

	case 2:
		*(Uint16*)p = pixel; // Set pixel for 16-bit surface
		break;

	case 3:
		if (SDL_BYTEORDER == SDL_BIG_ENDIAN) {
			p[0] = (pixel >> 16) & 0xff; // Set pixel for 24-bit surface on big-endian machine
			p[1] = (pixel >> 8) & 0xff;
			p[2] = pixel & 0xff;
		}
		else {
			p[0] = pixel & 0xff; // Set pixel for 24-bit surface on little-endian machine
			p[1] = (pixel >> 8) & 0xff;
			p[2] = (pixel >> 16) & 0xff;
		}
		break;

	case 4:
		*(Uint32*)p = pixel; // Set pixel for 32-bit surface
		break;
	}
}

// method to ensure colours don’t go out of 8 bit range in RGB​
void clamp255(Vec3f& col) {
	col.x = (col.x > 255) ? 255 : (col.x < 0) ? 0 : col.x; // Check if x is greater than 255 or less than 0 and clamp it to the range
	col.y = (col.y > 255) ? 255 : (col.y < 0) ? 0 : col.y; // Check if y is greater than 255 or less than 0 and clamp it to the range
	col.z = (col.z > 255) ? 255 : (col.z < 0) ? 0 : col.z; // Check if z is greater than 255 or less than 0 and clamp it to the range
}

// Finds the t-value of the intersection of a ray with a sphere.
double hit_sphere(const Point3f& centre, double radius, const Ray& r) {
	Vec3f oc = r.origin() - centre;
	auto a = r.direction().dotProduct(r.direction()); // Calculate dot product of ray direction vector
	auto b = 2.0 * oc.dotProduct(r.direction()); // Calculate dot product of the vector between ray origin and sphere center and ray direction vector
	auto c = oc.dotProduct(oc) - radius * radius; // Calculate the squared length of the vector between the ray origin and sphere center, subtracting the squared radius of the sphere
	auto discriminant = b * b - 4 * a * c; // Calculate the discriminant of the quadratic equation
	if (discriminant < 0) { // If discriminant is less than zero, the ray missed the sphere
		return -1.0;
	}
	else { // Otherwise, return the smaller t-value of the two intersections
		return (-b - sqrt(discriminant)) / (2.0 * a);
	}
}

// Determines the colour seen along a ray using recursive ray tracing.
Colour ray_colour(const Ray& r, const Colour& background, const hittable& world, int depth) {
	hit_record rec;
	// If we've exceeded the ray bounce limit, no more light is gathered.
	if (depth <= 0) return Colour(0, 0, 0); // Check if depth is less than or equal to zero, return black

	// If the ray hits nothing, return the background color.
	if (!world.hit(r, 0.001, infinity, rec))
		return background;

	Ray scattered;
	Colour attenuation;
	Colour emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

	if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
		return emitted;

	return emitted + attenuation * ray_colour(scattered, background, world, depth - 1);
}

// It takes in a pointer to the SDL_Surface that represents the screen,
// a hittable_list of objects in the world, the y-coordinate of the line to render,
// the number of samples per pixel, the maximum depth of recursion for rays, and a camera pointer
void LineRender(SDL_Surface* screen, hittable_list world, int y, int spp, int max_depth, camera* cam) {
	const float aspect_ratio = 16.0 / 9; // Set the aspect ratio of the image
	const int image_width = screen->w; // Get the width of the screen
	const int image_height = static_cast<int>(image_width / aspect_ratio); // Calculate the height of the image based on the aspect ratio

	const Colour black(0, 0, 0); // Create a black Colour object
	Colour pix_col(black); // Create a new Colour object and set it to black

	Colour background(0, 0, 0); // Create a black background Colour object

	// Loop through each pixel in the line
	for (int x = 0; x < screen->w; ++x) {
		pix_col = black; // Set the pixel Colour object to black
		// Loop through each sample for the pixel
		for (int s = 0; s < spp; s++) {
			// Calculate the u and v coordinates for the ray
			auto u = double(x + random_double()) / (image_width - 1);
			auto v = double(y + random_double()) / (image_height - 1);
			// Get the ray using the camera
			Ray ray = cam->get_ray(u, v);
			// Get the unit direction of the ray
			Vec3f unit_direction = ray.direction().normalize();
			// Calculate the background colour for the current pixel
			auto t = 0.5 * (unit_direction.y + 1.0);

			//background = (1.0 - t) * Colour(1.0, 1.0, 1.0) + t * Colour(135.0 / 255, 206.0 / 255, 250.0 / 255) * 255;
			//background = (1.0 - t) * Colour(1.0, 1.0, 1.0) + t * Colour(54.0 / 255, 142.0 / 255, 60.0 / 255) * 255; //background gradient for black into yellowish-turquoise
			//background = (1.0 - t) * Colour(1.0, 1.0, 1.0) + t * Colour(4.0 / 255, 92.0 / 255, 90.0 / 255) * 255; //background gradient for black into turquoise 
			// Get the colour of the pixel from tracing the ray
			pix_col = pix_col + ray_colour(ray, background, world, max_depth);
		}
		// Divide the pixel colour by the number of samples and normalize its brightness
		pix_col /= 255.f * spp;
		pix_col.x = sqrt(pix_col.x);
		pix_col.y = sqrt(pix_col.y);
		pix_col.z = sqrt(pix_col.z);
		pix_col *= 255;
		Uint32 colour = SDL_MapRGB(screen->format, pix_col.x, pix_col.y, pix_col.z);
		putpixel(screen, x, y, colour);
	}
}

//trying to render an earth using a texture
hittable_list earth() {
	auto earth_texture = make_shared<image_texture>("./textures/earthmap.jpg");
	auto earth_surface = make_shared<lambertian>(earth_texture);
	auto globe = make_shared<Sphere>(Point3f(0, 0, 0), 2, earth_surface);

	return hittable_list(globe);
}

hittable_list low_poly_tri_test() {
	hittable_list world;

	Model* model = new Model("./model/jelly.obj");
	//Model* model = new Model("./model/cc_t.obj");
	Model* model3 = new Model("./model/fish.obj");

	// Define a translation vector
	Vec3f transform(0, 0, 0);
	// Define a glass material with refraction index of 1.5
	auto jelly_light = make_shared<diffuse_light>(Colour(255, 243, 200));
	//auto jelly_light = make_shared<diffuse_light>(Colour(255 / 255.0, 243 / 255.0, 200 / 255.0));
	//auto material1 = make_shared<lambertian>(Colour(0.4, 0.8, 0.4));
	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model->vert(model->face(i)[0]);
		const Vec3f& v1 = model->vert(model->face(i)[1]);
		const Vec3f& v2 = model->vert(model->face(i)[2]);

		const Vec3f v0n = model->vn(model->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model->vn(model->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model->vn(model->vertex_normal(i)[2]);  //vertex normal for v2

		// Create a triangle with the vertices, transformed by the translation vector, and with the glass material
		world.add(make_shared<triangle>(v0 + transform, v1 + transform, v2 + transform, v0n, v1n, v2n, jelly_light));
	}

	// Update the translation vector for the next set of triangles
	Vec3f transform3 = Vec3f(-1.2, 0.8, 0);
	// Define a metal material with a gray color and 0.5 roughness
	auto mat_metal = make_shared<metal>(Colour(0.7, 0.6, 0.5), 0.5);
	// define the diffuse material with a color
	auto mat_diffuse = make_shared<lambertian>(Colour(0.8, 0.4, 0.7));
	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model3->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model3->vert(model3->face(i)[0]);
		const Vec3f& v1 = model3->vert(model3->face(i)[1]);
		const Vec3f& v2 = model3->vert(model3->face(i)[2]);

		const Vec3f v0n = model3->vn(model3->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model3->vn(model3->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model3->vn(model3->vertex_normal(i)[2]);  //vertex normal for v2

		world.add(make_shared<triangle>(v0 + transform3, v1 + transform3, v2 + transform3, v0n, v1n, v2n, mat_metal));
	}

	return hittable_list(make_shared<bvh_node>(world)); // with bvh
}

// This function generates a random scene by adding different types of spheres to a world object
hittable_list random_scene() {
	hittable_list world;

	// Add a ground sphere with a lambertian material to the world
	auto ground_material = make_shared<lambertian>(Colour(0.5, 0.5, 0.5));
	world.add(make_shared<Sphere>(Point3f(0, -1000, 0), 1000, ground_material));

	//Create a checker pattern sphere and add it to the world (currently commented out)
	auto checker = make_shared<checker_texture>(Colour(0.2, 0.3, 0.1), Colour(0.9, 0.9, 0.9));
	world.add(make_shared<Sphere>(Point3f(0, -1000, 0), 1000, make_shared<lambertian>(checker)));

	// Add multiple spheres with different materials to the world
	for (int a = -11; a < 11; a++) {
		for (int b = -11; b < 11; b++) {
			auto choose_mat = random_double();
			Point3f centre(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());
			if ((centre - Point3f(4, 0.2, 0)).length() > 0.9) {
				shared_ptr<material> sphere_material;

				if (choose_mat < 0.8) { // Diffuse material
					auto albedo = Colour::random() * Colour::random();
					sphere_material = make_shared<lambertian>(albedo);
					auto centre2 = centre + Vec3f(0, random_double(0, .5), 0);
					//add a moving sphere 
					world.add(make_shared<moving_sphere>(
						centre, centre2, 0.0, 1.0, 0.2, sphere_material));
					// Add a stationary sphere to the world (currently commented out)
					// world.add(make_shared<Sphere>(centre, 0.2, sphere_material));
				}
				else if (choose_mat < 0.95) { // Metal material
					auto albedo = Colour::random(0.5, 1);
					auto fuzz = random_double(0, 0.5);
					sphere_material = make_shared<metal>(albedo, fuzz);
					world.add(make_shared<Sphere>(centre, 0.2, sphere_material));
				}
				else { // Glass material
					sphere_material = make_shared<dielectric>(1.5);
					world.add(make_shared<Sphere>(centre, 0.2, sphere_material));
				}
			}
		}
	}

	// Add three more spheres with different materials to the world
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

// Define a hittable list called "world"
hittable_list final_scene() {
	hittable_list world;

	// Load a 3D model file into a "Model" object
	// Choose one of the following models to load by commenting/uncommenting the appropriate line
	//Model* model = new Model("./model/cc_t.obj");
	Model* model = new Model("./model/jelly.obj");
	Model* model2 = new Model("./model/coral.obj");
	Model* model3 = new Model("./model/fish.obj");
	Model* model4 = new Model("./model/foliage.obj");
	Model* model5 = new Model("./model/background.obj");

	// Define a translation vector
	Vec3f transform(0, 0, 0);
	// Define a glass material with refraction index of 1.5
	auto jelly_light = make_shared<diffuse_light>(Colour(255, 220, 180));
	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model->vert(model->face(i)[0]);
		const Vec3f& v1 = model->vert(model->face(i)[1]);
		const Vec3f& v2 = model->vert(model->face(i)[2]);

		const Vec3f v0n = model->vn(model->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model->vn(model->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model->vn(model->vertex_normal(i)[2]);  //vertex normal for v2

		// Create a triangle with the vertices, transformed by the translation vector, and with the glass material
		world.add(make_shared<triangle>(v0 + transform, v1 + transform, v2 + transform, v0n, v1n, v2n, jelly_light));
	}


	// Define a translation vector
	Vec3f transform2(0, 0, 0);
	shared_ptr<texture> coral_texture;
	coral_texture = make_shared<image_texture>("./textures/coral.jpg");
	auto coral_surface = make_shared<lambertian>(coral_texture);
	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model2->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model2->vert(model2->face(i)[0]);
		const Vec3f& v1 = model2->vert(model2->face(i)[1]);
		const Vec3f& v2 = model2->vert(model2->face(i)[2]);

		const Vec3f v0n = model2->vn(model2->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model2->vn(model2->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model2->vn(model2->vertex_normal(i)[2]);  //vertex normal for v2

		// Create a triangle with the vertices, transformed by the translation vector, and with the coral texture
		world.add(make_shared<triangle>(v0 + transform2, v1 + transform2, v2 + transform2, v0n, v1n, v2n, coral_surface));
	}


	// Update the translation vector for the next set of triangles
	Vec3f transform3 = Vec3f(-1.2, 0.8, 0);
	// Define a metal material with a gray color and 0.5 roughness
	auto mat_metal = make_shared<metal>(Colour(0.7, 0.6, 0.5), 0.5);
	// define the diffuse material with a color
	auto mat_diffuse = make_shared<lambertian>(Colour(0.8, 0.4, 0.7));
	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model3->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model3->vert(model3->face(i)[0]);
		const Vec3f& v1 = model3->vert(model3->face(i)[1]);
		const Vec3f& v2 = model3->vert(model3->face(i)[2]);

		const Vec3f v0n = model3->vn(model3->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model3->vn(model3->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model3->vn(model3->vertex_normal(i)[2]);  //vertex normal for v2

		world.add(make_shared<triangle>(v0 + transform3, v1 + transform3, v2 + transform3, v0n, v1n, v2n, mat_metal));
	}

	// Define a translation vector
	Vec3f transform4(0, 0, 0);
	// Define a diffuse material with a green color
	auto mat_diffuse2 = make_shared<lambertian>(Colour(0.1, 0.8, 0.1));
	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model4->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model4->vert(model4->face(i)[0]);
		const Vec3f& v1 = model4->vert(model4->face(i)[1]);
		const Vec3f& v2 = model4->vert(model4->face(i)[2]);

		const Vec3f v0n = model4->vn(model4->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model4->vn(model4->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model4->vn(model4->vertex_normal(i)[2]);  //vertex normal for v2

		// Create a triangle with the vertices, transformed by the translation vector, and with the diffuse material
		world.add(make_shared<triangle>(v0 + transform, v1 + transform, v2 + transform, v0n, v1n, v2n, mat_diffuse2));
	}

	// Define a translation vector
	Vec3f transform5(0, 0, 0);
	shared_ptr<texture> background_texture;
	background_texture = make_shared<image_texture>("./textures/underwater.jpg");
	auto background_surface = make_shared<lambertian>(background_texture);
	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model5->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model5->vert(model5->face(i)[0]);
		const Vec3f& v1 = model5->vert(model5->face(i)[1]);
		const Vec3f& v2 = model5->vert(model5->face(i)[2]);

		const Vec3f v0n = model5->vn(model5->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model5->vn(model5->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model5->vn(model5->vertex_normal(i)[2]);  //vertex normal for v2

		// Create a triangle with the vertices, transformed by the translation vector, and with texture
		world.add(make_shared<triangle>(v0 + transform5, v1 + transform5, v2 + transform5, v0n, v1n, v2n, background_surface));
	}


	// Add multiple spheres with different materials to the world
	for (int i = 0; i < 70; i++) {
		auto choose_mat = random_double();
		Point3f center(random_double(-20, 20), random_double(-0, 40), random_double(-20, 20));
		if ((center - Point3f(4, 0.2, 0)).length() > 0.9) {
			shared_ptr<material> sphere_material;

			if (choose_mat < 0.6) { // Diffuse lighting material
				sphere_material = make_shared<diffuse_light>(Colour(255, 243, 200)); // Warm yellow color
				world.add(make_shared<moving_sphere>(center, center + Vec3f(0, random_double(0, 2), 0), 0.0, 1.0, 0.5, sphere_material));
			}
			else if (choose_mat < 0.95 && (center - Point3f(2, 0.2, 0)).length() > 0.9) { // Metal material
				auto albedo = Colour(255 / 255.0, 243 / 255.0, 200 / 255.0);; // Warm yellow color
				auto fuzz = random_double(0, 0.5);
				sphere_material = make_shared<metal>(albedo, fuzz);
				world.add(make_shared<Sphere>(center, 1.0, sphere_material));
			}
			else { // Glass material
				sphere_material = make_shared<dielectric>(1.5);
				world.add(make_shared<Sphere>(center, 10.0, sphere_material));
			}
		}
	}

	//sandy ground
	auto ground_material = make_shared<lambertian>(Colour(194 / 255.0, 178 / 255.0, 128 / 255.0));
	world.add(make_shared<Sphere>(Point3f(0, -1000, 0), 1000, ground_material));
	auto light = make_shared<diffuse_light>(Colour(255, 255, 255));
	world.add(make_shared<Sphere>(Point3f(0, 100, 0), 10.0, light));

	//return world;                                       // without bvh
	return hittable_list(make_shared<bvh_node>(world)); // with bvh
}

hittable_list two_spheres() {
	hittable_list objects;

	auto checker = make_shared<checker_texture>(Colour(0.2, 0.3, 0.1), Colour(0.9, 0.9, 0.9));

	objects.add(make_shared<Sphere>(Point3f(0, -10, 0), 10, make_shared<lambertian>(checker)));
	objects.add(make_shared<Sphere>(Point3f(0, 10, 0), 10, make_shared<lambertian>(checker)));

	return objects;
}

hittable_list two_perlin_spheres() {
	hittable_list objects;
	auto pertext = make_shared<noise_texture>(4);
	objects.add(make_shared<Sphere>(Point3f(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
	objects.add(make_shared<Sphere>(Point3f(0, 2, 0), 2, make_shared<lambertian>(pertext)));

	return objects; //without bvh
	return hittable_list(make_shared<bvh_node>(objects)); // with bvh
}

hittable_list simple_light() {
	hittable_list objects;

	auto pertext = make_shared<noise_texture>(4);
	objects.add(make_shared<Sphere>(Point3f(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
	objects.add(make_shared<Sphere>(Point3f(0, 2, 0), 2, make_shared<lambertian>(pertext)));

	auto light = make_shared<diffuse_light>(Colour(255, 255, 255));
	objects.add(make_shared<Sphere>(Point3f(0, 100, 0), 10.0, light));

	auto difflight = make_shared<diffuse_light>(Colour(255, 255, 255));
	objects.add(make_shared<xy_rect>(3, 5, 1, 3, -2, difflight));

	return objects;
}

hittable_list cornell_box() {
	hittable_list objects;

	auto red = make_shared<lambertian>(Colour(.65, .05, .05));
	auto white = make_shared<lambertian>(Colour(.73, .73, .73));
	auto green = make_shared<lambertian>(Colour(.12, .45, .15));
	auto light = make_shared<diffuse_light>(Colour(255, 255, 255));

	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
	objects.add(make_shared<xz_rect>(213, 343, 227, 332, 554, light));
	objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
	objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
	objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

	//add instances for boxes
	shared_ptr<hittable> box1 = make_shared<box>(Point3f(0, 0, 0), Point3f(165, 330, 165), white);
	box1 = make_shared<rotate_y>(box1, 15);
	box1 = make_shared<translate>(box1, Vec3f(265, 0, 295));
	objects.add(box1);

	shared_ptr<hittable> box2 = make_shared<box>(Point3f(0, 0, 0), Point3f(165, 165, 165), white);
	box2 = make_shared<rotate_y>(box2, -18);
	box2 = make_shared<translate>(box2, Vec3f(130, 0, 65));
	objects.add(box2);
	return objects;
}

hittable_list cornell_smoke() {
	hittable_list objects;

	auto red = make_shared<lambertian>(Colour(.65, .05, .05));
	auto white = make_shared<lambertian>(Colour(.73, .73, .73));
	auto green = make_shared<lambertian>(Colour(.12, .45, .15));
	auto light = make_shared<diffuse_light>(Colour(7, 7, 7));

	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
	objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
	objects.add(make_shared<xz_rect>(113, 443, 127, 432, 554, light));
	objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
	objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
	objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

	shared_ptr<hittable> box1 = make_shared<box>(Point3f(0, 0, 0), Point3f(165, 330, 165), white);
	box1 = make_shared<rotate_y>(box1, 15);
	box1 = make_shared<translate>(box1, Vec3f(265, 0, 295));

	shared_ptr<hittable> box2 = make_shared<box>(Point3f(0, 0, 0), Point3f(165, 165, 165), white);
	box2 = make_shared<rotate_y>(box2, -18);
	box2 = make_shared<translate>(box2, Vec3f(130, 0, 65));

	objects.add(make_shared<constant_medium>(box1, 0.01, Colour(0, 0, 0)));
	objects.add(make_shared<constant_medium>(box2, 0.01, Colour(1, 1, 1)));

	return objects;
}


hittable_list combined_scene() {
	hittable_list boxes1;
	auto ground = make_shared<lambertian>(Colour(0.48, 0.83, 0.53));

	const int boxes_per_side = 20;
	for (int i = 0; i < boxes_per_side; i++) {
		for (int j = 0; j < boxes_per_side; j++) {
			auto w = 100.0;
			auto x0 = -1000.0 + i * w;
			auto z0 = -1000.0 + j * w;
			auto y0 = 0.0;
			auto x1 = x0 + w;
			auto y1 = random_double(1, 101);
			auto z1 = z0 + w;

			boxes1.add(make_shared<box>(Point3f(x0, y0, z0), Point3f(x1, y1, z1), ground));
		}
	}

	hittable_list objects;

	objects.add(make_shared<bvh_node>(boxes1));

	auto light = make_shared<diffuse_light>(Colour(7, 7, 7));
	objects.add(make_shared<xz_rect>(123, 423, 147, 412, 554, light));

	auto center1 = Point3f(400, 400, 200);
	auto center2 = center1 + Vec3f(30, 0, 0);
	auto moving_sphere_material = make_shared<lambertian>(Colour(0.7, 0.3, 0.1));
	objects.add(make_shared<moving_sphere>(center1, center2, 0, 1, 50, moving_sphere_material));

	objects.add(make_shared<Sphere>(Point3f(260, 150, 45), 50, make_shared<dielectric>(1.5)));
	objects.add(make_shared<Sphere>(
		Point3f(0, 150, 145), 50, make_shared<metal>(Colour(0.8, 0.8, 0.9), 1.0)
		));

	auto boundary = make_shared<Sphere>(Point3f(360, 150, 145), 70, make_shared<dielectric>(1.5));
	objects.add(boundary);
	objects.add(make_shared<constant_medium>(boundary, 0.2, Colour(0.2, 0.4, 0.9)));
	boundary = make_shared<Sphere>(Point3f(0, 0, 0), 5000, make_shared<dielectric>(1.5));
	objects.add(make_shared<constant_medium>(boundary, .0001, Colour(1, 1, 1)));

	auto emat = make_shared<lambertian>(make_shared<image_texture>("./textures/earthmap.jpg"));
	objects.add(make_shared<Sphere>(Point3f(400, 200, 400), 100, emat));
	auto pertext = make_shared<noise_texture>(0.1);
	objects.add(make_shared<Sphere>(Point3f(220, 280, 300), 80, make_shared<lambertian>(pertext)));

	hittable_list boxes2;
	auto white = make_shared<lambertian>(Colour(.73, .73, .73));
	int ns = 1000;
	for (int j = 0; j < ns; j++) {
		boxes2.add(make_shared<Sphere>(Point3f::random(0, 165), 10, white));
	}

	objects.add(make_shared<translate>(
		make_shared<rotate_y>(
			make_shared<bvh_node>(boxes2), 15),
		Vec3f(-100, 270, 395)
		)
	);

	return objects;
}

hittable_list space_scene() {
	hittable_list world;

	// Load a 3D model file into a "Model" object
	// Choose one of the following models to load by commenting/uncommenting the appropriate line
	//Model* model = new Model("./model/cc_t.obj");
	Model* model = new Model("./model/astronaunt.obj");
	Model* model2 = new Model("./model/ground2.obj");
	Model* model3 = new Model("./model/asteroid.obj");
	Model* model4 = new Model("./model/planet.obj");
	Model* model5 = new Model("./model/backgroundSpace.obj");
	Model* model6 = new Model("./model/blackhole_inside.obj");
	Model* model7 = new Model("./model/blackhole_outside1.obj");
	Model* model8 = new Model("./model/blackhole_outside2.obj");
	Model* model9 = new Model("./model/blackhole_outside3.obj");
	Model* model10 = new Model("./model/blackhole_ring.obj");

	// Define a translation vector
	Vec3f transform(0, 0, 0);
	// Define a diffuse material with a white color
	//auto mat_diffuse = make_shared<lambertian>(Colour(255, 255, 255));
	auto AS_light = make_shared<diffuse_light>(Colour(255, 255, 255));
	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model->vert(model->face(i)[0]);
		const Vec3f& v1 = model->vert(model->face(i)[1]);
		const Vec3f& v2 = model->vert(model->face(i)[2]);

		const Vec3f v0n = model->vn(model->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model->vn(model->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model->vn(model->vertex_normal(i)[2]);  //vertex normal for v2

		// Create a triangle with the vertices, transformed by the translation vector, and with the diffuse material
		world.add(make_shared<triangle>(v0 + transform, v1 + transform, v2 + transform, v0n, v1n, v2n, AS_light));
	}

	// Update the translation vector for the next set of triangles
	Vec3f transform2 = Vec3f(-1.2, 0.8, 0);
	// Define a metal material with a gray color and 0.2 roughness
	auto mat_metal = make_shared<metal>(Colour(0.7, 0.6, 0.5), 0.2);
	//auto mat_metal = make_shared<metal>(Colour(0.7, 0.6, 0.5), 0.2);
	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model2->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model2->vert(model2->face(i)[0]);
		const Vec3f& v1 = model2->vert(model2->face(i)[1]);
		const Vec3f& v2 = model2->vert(model2->face(i)[2]);

		const Vec3f v0n = model2->vn(model2->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model2->vn(model2->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model2->vn(model2->vertex_normal(i)[2]);  //vertex normal for v2

		world.add(make_shared<triangle>(v0 + transform2, v1 + transform2, v2 + transform2, v0n, v1n, v2n, mat_metal));
	}

	// Update the translation vector for the next set of triangles
	Vec3f transform3 = Vec3f(-1.2, 0.8, 0);
	// Define a metal material with a gray color and 0.5 roughness
	auto mat_metal2 = make_shared<metal>(Colour(0.7, 0.6, 0.5), 0.5);
	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model3->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model3->vert(model3->face(i)[0]);
		const Vec3f& v1 = model3->vert(model3->face(i)[1]);
		const Vec3f& v2 = model3->vert(model3->face(i)[2]);

		const Vec3f v0n = model3->vn(model3->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model3->vn(model3->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model3->vn(model3->vertex_normal(i)[2]);  //vertex normal for v2

		world.add(make_shared<triangle>(v0 + transform3, v1 + transform3, v2 + transform3, v0n, v1n, v2n, mat_metal2));
	}

	// Define a translation vector
	Vec3f transform4(0, 0, 0);
	shared_ptr<texture> planet_texture;
	planet_texture = make_shared<image_texture>("./textures/planet.jpg");
	auto planet_surface = make_shared<lambertian>(planet_texture);
	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model4->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model4->vert(model4->face(i)[0]);
		const Vec3f& v1 = model4->vert(model4->face(i)[1]);
		const Vec3f& v2 = model4->vert(model4->face(i)[2]);

		const Vec3f v0n = model4->vn(model4->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model4->vn(model4->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model4->vn(model4->vertex_normal(i)[2]);  //vertex normal for v2

		// Create a triangle with the vertices, transformed by the translation vector, and with the texture
		world.add(make_shared<triangle>(v0 + transform4, v1 + transform4, v2 + transform4, v0n, v1n, v2n, planet_surface));
	}

	// Define a diffuse material with a black color
	auto mat_diffuse3 = make_shared<lambertian>(Colour(0, 0, 0));

	//// Define a translation vector
	//Vec3f transform5(0, 0, 0);
	//// Loop through each triangle face in the loaded model
	//for (uint32_t i = 0; i < model5->nfaces(); ++i) {
	//	// Extract the three vertices of the triangle face
	//	const Vec3f& v0 = model5->vert(model5->face(i)[0]);
	//	const Vec3f& v1 = model5->vert(model5->face(i)[1]);
	//	const Vec3f& v2 = model5->vert(model5->face(i)[2]);

	//	const Vec3f v0n = model5->vn(model5->vertex_normal(i)[0]);  //vertex normal for v0
	//	const Vec3f v1n = model5->vn(model5->vertex_normal(i)[1]);  //vertex normal for v1
	//	const Vec3f v2n = model5->vn(model5->vertex_normal(i)[2]);  //vertex normal for v2

	//	// Create a triangle with the vertices, transformed by the translation vector, and with the texture
	//	world.add(make_shared<triangle>(v0 + transform5, v1 + transform5, v2 + transform5, v0n, v1n, v2n, mat_diffuse3));
	//}
	

	auto BH_light = make_shared<diffuse_light>(Colour(255, 255, 255));

	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model6->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model6->vert(model6->face(i)[0]);
		const Vec3f& v1 = model6->vert(model6->face(i)[1]);
		const Vec3f& v2 = model6->vert(model6->face(i)[2]);

		const Vec3f v0n = model6->vn(model6->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model6->vn(model6->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model6->vn(model6->vertex_normal(i)[2]);  //vertex normal for v2

		// Create a triangle with the vertices, transformed by the translation vector, and with the texture
		world.add(make_shared<triangle>(v0, v1, v2, v0n, v1n, v2n, mat_diffuse3));
	}	
	
	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model7->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model7->vert(model7->face(i)[0]);
		const Vec3f& v1 = model7->vert(model7->face(i)[1]);
		const Vec3f& v2 = model7->vert(model7->face(i)[2]);

		const Vec3f v0n = model7->vn(model7->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model7->vn(model7->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model7->vn(model7->vertex_normal(i)[2]);  //vertex normal for v2

		// Create a triangle with the vertices, transformed by the translation vector, and with the texture
		world.add(make_shared<triangle>(v0, v1, v2, v0n, v1n, v2n, BH_light));
	}

	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model8->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model8->vert(model8->face(i)[0]);
		const Vec3f& v1 = model8->vert(model8->face(i)[1]);
		const Vec3f& v2 = model8->vert(model8->face(i)[2]);

		const Vec3f v0n = model8->vn(model8->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model8->vn(model8->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model8->vn(model8->vertex_normal(i)[2]);  //vertex normal for v2

		// Create a triangle with the vertices, transformed by the translation vector, and with the coral texture
		world.add(make_shared<triangle>(v0, v1, v2, v0n, v1n, v2n, mat_diffuse3));
	}

	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model9->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model9->vert(model9->face(i)[0]);
		const Vec3f& v1 = model9->vert(model9->face(i)[1]);
		const Vec3f& v2 = model9->vert(model9->face(i)[2]);

		const Vec3f v0n = model9->vn(model9->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model9->vn(model9->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model9->vn(model9->vertex_normal(i)[2]);  //vertex normal for v2

		// Create a triangle with the vertices, transformed by the translation vector, and with the coral texture
		world.add(make_shared<triangle>(v0, v1, v2, v0n, v1n, v2n, BH_light));
	}

	// Loop through each triangle face in the loaded model
	for (uint32_t i = 0; i < model10->nfaces(); ++i) {
		// Extract the three vertices of the triangle face
		const Vec3f& v0 = model10->vert(model10->face(i)[0]);
		const Vec3f& v1 = model10->vert(model10->face(i)[1]);
		const Vec3f& v2 = model10->vert(model10->face(i)[2]);

		const Vec3f v0n = model10->vn(model10->vertex_normal(i)[0]);  //vertex normal for v0
		const Vec3f v1n = model10->vn(model10->vertex_normal(i)[1]);  //vertex normal for v1
		const Vec3f v2n = model10->vn(model10->vertex_normal(i)[2]);  //vertex normal for v2

		// Create a triangle with the vertices, transformed by the translation vector, and with the coral texture
		world.add(make_shared<triangle>(v0, v1, v2, v0n, v1n, v2n, BH_light));
	}

	// Generate spheres around lookat point
	for (int i = 0; i < 1000; i++) {
		auto choose_mat = random_double();
		Point3f center = lookat + Point3f(random_double(-1000, 1000), random_double(-700, 700), random_double(-1300, -1000));
		if ((center - lookat).length() > 0.9) {
			shared_ptr<material> sphere_material;

			if (choose_mat < 0.6) { // Diffuse lighting material
				auto albedo = Colour::random() * Colour::random();
				sphere_material = make_shared<lambertian>(albedo);
				world.add(make_shared<moving_sphere>(center, center + Vec3f(0, random_double(0, 2), 0), 0.0, 1.0, 2.f, sphere_material));
			}
			else if (choose_mat < 0.95 && (center - Point3f(2, 0.2, 0)).length() > 0.9) { // Metal material
				auto albedo = Colour::random() * Colour::random();
				auto fuzz = random_double(0, 0.5);
				sphere_material = make_shared<metal>(albedo, fuzz);
				world.add(make_shared<Sphere>(center, 2.f, sphere_material));
			}
			else { // Glass material
				sphere_material = make_shared<dielectric>(1.5);
				world.add(make_shared<Sphere>(center, 2.f, sphere_material));
			}
		}
	}

	auto purple = make_shared<lambertian>(Colour(.6, .1, .6));

	// Scale factor for the boxes
	float scaleXYZ = 10.0f;  // Significantly increase the size in all dimensions

	// Random ranges for box positions in the X and Y axes
	//float randX = random_double(-1000, 1000);
	float randX = -400;
	float fixedY = 0;
	float fixedZ = -4300;  

	/*shared_ptr<hittable> box1 = make_shared<box>(Point3f(0, 0, 0), Point3f(165 * scaleXYZ, 330 * scaleXYZ, 165 * scaleXYZ), purple);
	box1 = make_shared<rotate_y>(box1, 15);
	box1 = make_shared<translate>(box1, lookat + Vec3f(randX, fixedY, fixedZ));*/

	shared_ptr<hittable> box2 = make_shared<box>(Point3f(0, 0, 0), Point3f(300 * scaleXYZ, 100 * scaleXYZ, 100), purple);
	box2 = make_shared<rotate_y>(box2, -18);
	box2 = make_shared<translate>(box2, lookat + Vec3f(randX, fixedY, fixedZ));

	//world.add(make_shared<constant_medium>(box1, 0.01, Colour(0, 0, 0)));
	world.add(make_shared<constant_medium>(box2, 1.f, Colour(.6, .1, .6)));

	//return world;                                       // without bvh
	return hittable_list(make_shared<bvh_node>(world)); // with bvh
}


Uint32 getpixel(SDL_Surface* surface, int x, int y)
{
	int bpp = surface->format->BytesPerPixel;
	/* Here p is the address to the pixel we want to retrieve */
	Uint8* p = (Uint8*)surface->pixels + y * surface->pitch + x * bpp;

	switch (bpp)
	{
	case 1:
		return *p;
		break;

	case 2:
		return *(Uint16*)p;
		break;

	case 3:
		if (SDL_BYTEORDER == SDL_BIG_ENDIAN)
			return p[0] << 16 | p[1] << 8 | p[2];
		else
			return p[0] | p[1] << 8 | p[2] << 16;
		break;

	case 4:
		return *(Uint32*)p;
		break;

	default:
		return 0;       /* shouldn't happen, but avoids warnings */
	}
}

//initialise
int main(int argc, char** argv)
{
	init();

	playMusic();

	// Image settings
	int spp = 1; // samples per pixel
	const auto aspect_ratio = 16.0 / 9.0; // width to height ratio of the image
	const int image_width = screen->w; // width of the output image
	const int image_height = static_cast<int>(image_width / aspect_ratio); // height of the output image
	const float scale = 1.0f / spp; // scale the colour values according to the number of samples per pixel
	const int max_depth = 50; // maximum depth for recursive ray tracing

	// Create the scene
	switch (0) {

		// Default case for when no matching case is found
		// Set the camera position and view angle for each scene
		// The settings for each scene can be adjusted as needed
		//default:
	case 1:
		world = random_scene();
		lookfrom = Point3f(13, 2, 3);
		lookat = Point3f(0, 0, 0);
		vfov = 20.0;
		aperture = 0.1;
		vup = Vec3f(0, 10, 0);
		dist_to_focus = 10.0;
		background = Colour(0.70, 0.80, 1.00);
		//dist_to_focus = (lookfrom - lookat).length();
		break;

		//default:
	case 3:
		world = two_spheres();
		lookfrom = Point3f(13, 2, 3);
		lookat = Point3f(0, 0, 0);
		vfov = 20.0;
		background = Colour(0.70, 0.80, 1.00);
		break;

		//default: 
	case 2:
		world = final_scene();
		lookfrom = Point3f(-25.148f, 26.359f, 47.128f);
		lookat = Point3f(0.f, 17.f, 0.f);
		vfov = 35.000f;
		////aperture = 0.1;
		aperture = 0.7f;
		vup = Vec3f(0, 1, 0);
		background = Colour(0.70, 0.80, 1.00);
		dist_to_focus = (lookfrom - lookat).length();
		break;

		//default:
	case 4:
		world = earth();
		lookfrom = Point3f(13, 2, 3);
		lookat = Point3f(0, 0, 0);
		vfov = 20.0;
		break;

		//default:
	case 5:
		world = low_poly_tri_test();
		lookfrom = Point3f(-25.148f, 26.359f, 47.128f);
		lookat = Point3f(0.f, 17.f, 0.f);
		vfov = 35.000f;
		aperture = 0.7f;
		vup = Vec3f(0, 1, 0);
		background = Colour(0.70, 0.80, 1.00);
		dist_to_focus = (lookfrom - lookat).length();
		break;

		//default:
	case 6:
		world = two_perlin_spheres();
		vfov = 20.0;
		background = Colour(0.70, 0.80, 1.00);

		//for camera positions view the inner render loop

		break;

		//default:
	case 7:
		world = simple_light();
		spp = 1;
		lookfrom = Point3f(26, 3, 6);
		lookat = Point3f(0, 2, 0);
		background = Colour(0.0, 0.0, 0.0);
		vfov = 20.0;
		break;

		//default:
	case 8:
		background = Colour(0.0, 0.0, 0.0);
		break;

	//default:
	case 9:
		world = cornell_box();
		spp = 10;
		background = Colour(0, 0, 0);
		lookfrom = Point3f(278, 278, -800);
		lookat = Point3f(278, 278, 0);
		vfov = 40.0;
		break;

	//default:
	case 10:
		world = cornell_smoke();
		spp = 10;
		lookfrom = Point3f(278, 278, -800);
		lookat = Point3f(278, 278, 0);
		vfov = 40.0;
		break;

	//default:
	case 11:
		world = combined_scene();
		spp = 10;
		background = Colour(0, 0, 0);
		lookfrom = Point3f(478, 278, -600);
		lookat = Point3f(278, 278, 0);
		vfov = 40.0;
		break;

	default:
	case 12:
		world = space_scene();
		vfov = -1040.000f;
		//vfov = 35.000f;
		background = Colour(0, 0, 0);
		aperture = 1.417f;
		vup = Vec3f(0, 1, 0);
		spp = 1;
		lookfrom = Point3f(camX + -180.122f, camY + 133.769f, camZ + 202.863f);
		lookat = Point3f(-140, 150, 0);
		dist_to_focus = (lookfrom - lookat).length();
		
		//for camera positions view the inner render loop
		

		break;
	}

	// Initialize the camera with the specified parameters
	camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

	double t;
	Colour pix_col;
	Ray r;
	pix_col.x += ray_colour (r, background, world, max_depth).x;
	pix_col.y += ray_colour (r, background, world, max_depth).y;
	pix_col.z += ray_colour (r, background, world, max_depth).z;

	// Render loop
	SDL_Event e;
	bool running = true;
	while (running) {
		// Start a timer to measure the frame render time
		auto t_start = std::chrono::high_resolution_clock::now();

		// Calculate the new camera position based on the angle
		camAngleX += 0.01f;
		camAngleY += 0.01f;

		// Calculate the new camera position
		camX = distance * -sinf(camAngleX) * cosf(camAngleY);
		camY = distance * -sinf(camAngleY);
		camZ = distance * cosf(camAngleX) * cosf(camAngleY);

		// Update the lookfrom and lookat points
		lookfrom = Point3f(camX + -180.122f, camY + 133.769f, camZ + 202.863f);
		lookat = Point3f(-140, 150, 0);
		
		/*lookfrom = Point3f(camX + 13, camY + 2, camZ + 3);
		lookat = Point3f(0, 0, 0);*/

		// Update the camera in the rendering pipeline
		cam.updateCamera(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

		// Clear the back buffer, pixel data on the surface, and depth buffer
		// ...

		//multithreading
		{
			// Start a timer to measure the render time of each line
			t_start = std::chrono::high_resolution_clock::now();

			// Create a thread pool and divide the rendering workload
			ThreadPool pool(std::thread::hardware_concurrency());

			int start = screen->h - 1;
			int step = screen->h / std::thread::hardware_concurrency();

			// Enqueue each line of the screen to be rendered by the thread pool
			for (int y = 0; y < screen->h - 1; y++)
			{
				pool.Enqueue(std::bind(LineRender, screen, world, y, spp, max_depth, &cam));
			}
		}

		// Calculate the frame render time and print it out
		auto t_end = std::chrono::high_resolution_clock::now();
		auto passedTime = std::chrono::duration<double, std::milli>(t_end - t_start).count();
		std::cerr << "Frame render time:  " << passedTime << " ms" << std::endl;

		// Save the rendered image to a file in PPM format
		std::ofstream out("raycaster.ppm");
		out << "P3\n" << image_width << " " << image_height << "\n255\n";
		for (int j = image_height - 1; j >= 0; --j) {
			for (int i = 0; i < image_width; ++i) {
				SDL_Color rgb;
				Uint32 data = getpixel(screen, i, j);
				SDL_GetRGB(data, screen->format, &rgb.r, &rgb.g, &rgb.b);
				//std::cout << "rgb: " << static_cast<int>(rgb.r) << " " << static_cast<int>(rgb.g) << " " << static_cast<int>(rgb.b) << std::endl;
				out << static_cast<int>(rgb.r) << ' ' << static_cast<int>(rgb.g) << ' ' << static_cast<int>(rgb.b) << '\n';
			}
		}
		out.close();

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

		// Check for SDL events
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

	// Quit SDL_mixer
	Mix_Quit();

	return 0;
}