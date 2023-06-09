#include "geometry.h"
#include "SDL.h" 
#include "model.h"
#include "vertexShader.h"
#include "tgaimage.h"
#include <fstream>
#include <chrono>

#define M_PI 3.14159265359

struct PointLight {
    Vec3f position;
    Vec3f color;
    float intensity;
};

enum FitResolutionGate { kFill = 0, kOverscan };

//variables
const uint32_t imageWidth = 1920;
const uint32_t imageHeight = 1080;
Matrix44f worldToCamera;
static const float inchToMm = 25.4;
const float nearClippingPlane = 0.1;
const float farClippingPlane = 10000.f;
float focalLength = 35; // in mm
// 35mm Full Aperture in inches
float filmApertureWidth = 36.f;
float filmApertureHeight = 24.f;
// Initialize accumulators for color and depth values
Vec3f accumColour(0.f, 0.f, 0.f);
float accumDepth = 0;
int sampleCount = 0;
Vec3f lightDirection(0.f, -1.f, 0.f); // A light position in camera space
Vec3f lightColour(1.2f, 1.1f, 0.9f); // colours the light warm white
float shininess = 32.f; //shininess of the specular lighting reflection
TGAImage texture, texture2, texture3, texture4, texture5, texture6;
TGAImage specular, specular2, specular3, specular4, specular5, specular6;

//Pointers
SDL_Window* window;
SDL_Renderer* renderer;
SDL_Surface* screen;
Model* model = nullptr;
Model* model2 = nullptr;
Model* model3 = nullptr;
Model* model4 = nullptr;
Model* model5 = nullptr;
Model* model6 = nullptr;
Model* currentModel[6];

// Normalizes the given 3D vector in-place, so that its length becomes 1.
void normalize(Vec3f& v)
{
    // Compute the length of the vector using the Pythagorean theorem.
    float len = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    // Print the original and normalized vector components for debugging purposes.
    std::cout << v.x << " " << v.y << " " << v.z << "\n";
    // Divide each component of the vector by its length to obtain the normalized vector.
    v.x /= len, v.y /= len, v.z /= len;
    // Print the normalized vector components for debugging purposes.
    std::cout << v.x << " " << v.y << " " << v.z << "\n";
}

// Computes the cross product of two 3D vectors a and b, which is a vector
// that is perpendicular to both a and b, and has a magnitude equal to the
// area of the parallelogram formed by a and b.
Vec3f cross(const Vec3f& a, const Vec3f& b)
{
    // Compute the components of the cross product using the determinant formula.
    return {
    a.y * b.z - a.z * b.y, // x component
    a.z * b.x - a.x * b.z, // y component
    a.x * b.y - a.y * b.x // z component
    };
}


void computeScreenCoordinates(
    const float &filmApertureWidth,
    const float &filmApertureHeight,
    const uint32_t &imageWidth,
    const uint32_t &imageHeight,
    const FitResolutionGate &fitFilm,
    const float &nearClippingPlane,
    const float &focalLength,
    float &top, float &bottom, float &left, float &right
)
{
    float filmAspectRatio = filmApertureWidth / filmApertureHeight;
    float deviceAspectRatio = imageWidth / (float)imageHeight;
    
    top = ((filmApertureHeight * inchToMm / 2) / focalLength) * nearClippingPlane;
    right = ((filmApertureWidth * inchToMm / 2) / focalLength) * nearClippingPlane;

    // field of view (horizontal)
    //float fov = 2 * 180 / M_PI * atan((filmApertureWidth * inchToMm / 2) / focalLength);
    float fov = 54.4322231146;
    std::cerr << "Field of view " << fov << std::endl;
    
    float xscale = 1;
    float yscale = 1;
    
    switch (fitFilm) {
        default:
        case kFill:
            if (filmAspectRatio > deviceAspectRatio) {
                xscale = deviceAspectRatio / filmAspectRatio;
            }
            else {
                yscale = filmAspectRatio / deviceAspectRatio;
            }
            break;
        case kOverscan:
            if (filmAspectRatio > deviceAspectRatio) {
                yscale = filmAspectRatio / deviceAspectRatio;
            }
            else {
                xscale = deviceAspectRatio / filmAspectRatio;
            }
            break;
    }
    
    right *= xscale;
    top *= yscale;
    
    bottom = -top;
    left = -right;
}

// Compute vertex raster screen coordinates.
// Vertices are defined in world space. They are then converted to camera space,
// then to NDC space (in the range [-1,1]) and then to raster space.
// The z-coordinates of the vertex in raster space is set with the z-coordinate
// of the vertex in camera space.
void convertToRaster(
    const Vec3f &vertexWorld,
    const Matrix44f &worldToCamera,
    const float &l,
    const float &r,
    const float &t,
    const float &b,
    const float &near,
    const uint32_t &imageWidth,
    const uint32_t &imageHeight,
    Vec3f &vertexRaster
)
{
    // convert to camera space
    float nearClippingPlane = 1.;
    // point in camera space
    Vec3f vertexCamera;
    worldToCamera.multVecMatrix(vertexWorld, vertexCamera);
     
    // convert to screen space
    Vec2f vertexScreen;
    vertexScreen.x = nearClippingPlane * vertexCamera.x / -vertexCamera.z;
    vertexScreen.y = nearClippingPlane * vertexCamera.y / -vertexCamera.z;

    // now convert point from screen space to NDC space (in range [-1,1])
    Vec2f vertexNDC;
    vertexNDC.x = 2 * vertexScreen.x / (r - l) - (r + l) / (r - l);
    vertexNDC.y = 2 * vertexScreen.y / (t - b) - (t + b) / (t - b);
   
    // convert to raster space and set point z-coordinate to -vertexCamera.z
    vertexRaster.x = (vertexNDC.x + 1) / 2 * imageWidth;
    // in raster space y is down so invert direction
    vertexRaster.y = (1 - vertexNDC.y) / 2 * imageHeight;
    // store the point camera space z-coordinate (as a positive value)
    vertexRaster.z = -vertexCamera.z;
}

// Returns the minimum value among three given floating-point numbers.
float min3(const float& a, const float& b, const float& c)
{
    return std::min(a, std::min(b, c)); // using the standard library's min function
}

// Returns the maximum value among three given floating-point numbers.
float max3(const float& a, const float& b, const float& c)
{
    return std::max(a, std::max(b, c)); // using the standard library's max function
}

// Computes the edge function of three vertices a, b, and c, which is a determinant
// of a 2x2 matrix formed by subtracting a from c and b from a.
// The result is positive if c lies to the left of the line passing through a and b,
// negative if c lies to the right, and zero if c lies on the line.
float edgeFunction(const Vec3f& a, const Vec3f& b, const Vec3f& c)
{
    return (c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0]);
}

// Initialize SDL and create window and renderer
void init() {

    SDL_Init(SDL_INIT_VIDEO);

    // Create window
    window = SDL_CreateWindow(
        "Software Rasteriser",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        1920,
        1080,
        0
    );

    // Get window surface
    screen = SDL_GetWindowSurface(window);
    // Create renderer
    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
    // Set render colour to black
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
}


// Set the pixel for 8-bit per pixel surfaces.
void set_pixel_8bpp(Uint8* p, Uint32 pixel)
{
    *p = pixel;
}

// Set the pixel for 16-bit per pixel surfaces.
void set_pixel_16bpp(Uint8* p, Uint32 pixel)
{
    *(Uint16*)p = pixel;
}

// Set the pixel for 24-bit per pixel surfaces.
void set_pixel_24bpp(Uint8* p, Uint32 pixel)
{
    if (SDL_BYTEORDER == SDL_BIG_ENDIAN) {
        // Set the pixel for big-endian systems.
        p[0] = (pixel >> 16) & 0xff;
        p[1] = (pixel >> 8) & 0xff;
        p[2] = pixel & 0xff;
    }
    else {
        // Set the pixel for little-endian systems.
        p[0] = pixel & 0xff;
        p[1] = (pixel >> 8) & 0xff;
        p[2] = (pixel >> 16) & 0xff;
    }
}

// Set the pixel for 32-bit per pixel surfaces.
void set_pixel_32bpp(Uint8* p, Uint32 pixel)
{
    *(Uint32*)p = pixel;
}

// Set the pixel at position (x,y) on the surface to the specified pixel value.
void putpixel(SDL_Surface* surface, int x, int y, Uint32 pixel)
{
    int bpp = surface->format->BytesPerPixel;
    // Here p is the address to the pixel we want to set.
    Uint8* p = (Uint8*)surface->pixels + y * surface->pitch + x * bpp;

    switch (bpp) {
        // For 8 bits per pixel surfaces, set the pixel.
    case 1:
        set_pixel_8bpp(p, pixel);
        break;

        // For 16 bits per pixel surfaces, set the pixel.
    case 2:
        set_pixel_16bpp(p, pixel);
        break;

        // For 24 bits per pixel surfaces, set the pixel.
    case 3:
        set_pixel_24bpp(p, pixel);
        break;

        // For 32 bits per pixel surfaces, set the pixel.
    case 4:
        set_pixel_32bpp(p, pixel);
        break;
    }
}


// This function calculates the forward vector given two 3D vectors: from and to
Vec3f calculateForwardVector(const Vec3f from, const Vec3f to)
{
    // Calculate the forward vector
    Vec3f forward = from - to;
    // Normalize the forward vector
    normalize(forward);
    // Return the normalized forward vector
    return forward;
}

// This function calculates the right vector given two 3D vectors: up and forward
Vec3f calculateRightVector(const Vec3f up, const Vec3f forward)
{
    // Calculate the right vector using cross product
    Vec3f right = cross(up, forward);
    // Normalize the right vector
    normalize(right);
    // Return the normalized right vector
    return right;
}

// This function calculates the up vector given two 3D vectors: forward and right
Vec3f calculateUpVector(const Vec3f forward, const Vec3f right)
{
    // Calculate the up vector using cross product
    Vec3f up = cross(forward, right);
    // Return the up vector
    return up;
}

// This function generates a 4x4 transformation matrix that transforms from the camera space to world space
Matrix44f lookAt(const Vec3f from, const Vec3f to, const Vec3f up, const Vec3f _tmp = Vec3f(0, 1, 0))
{
    // Calculate the forward, right, and up vectors
    Vec3f forward = calculateForwardVector(from, to);
    Vec3f right = calculateRightVector(up, forward);
    Vec3f newup = calculateUpVector(forward, right);

    // Create the transformation matrix using the calculated vectors and the camera position
    Matrix44f camToWorld;
    camToWorld[0][0] = right.x, camToWorld[0][1] = right.y, camToWorld[0][2] = right.z;
    camToWorld[1][0] = newup.x, camToWorld[1][1] = newup.y, camToWorld[1][2] = newup.z;
    camToWorld[2][0] = forward.x, camToWorld[2][1] = forward.y, camToWorld[2][2] = forward.z;
    camToWorld[3][0] = from.x, camToWorld[3][1] = from.y, camToWorld[3][2] = from.z;

    // Return the transformation matrix
    return camToWorld;
}


int main(int argc, char **argv)
{
    if (2 == argc) { model = new Model(argv[1]); }
    else { 
        //model = new Model("./model/cc_t.obj"); 
        model = new Model("./model/jelly.obj"); 
        model2 = new Model("./model/foliage.obj");
        model3 = new Model("./model/ground_rast2.obj"); 
        model4 = new Model("./model/coral.obj"); 
        model5 = new Model("./model/fish.obj");
        model6 = new Model("./model/background.obj");
    }

    //textures for model
    texture.read_tga_file("./textures/jellytex.tga");
    texture.flip_vertically(); // Flip the texture vertically because TGA files are stored upside down
    
    texture2.read_tga_file("./textures/seaweed.tga");
    texture2.flip_vertically(); // Flip the texture vertically because TGA files are stored upside down

    texture3.read_tga_file("./textures/sand.tga");
    texture3.flip_vertically(); // Flip the texture vertically because TGA files are stored upside downTGAImage texture3;

    texture4.read_tga_file("./textures/coral.tga");
    texture4.flip_vertically(); // Flip the texture vertically because TGA files are stored upside down

    texture5.read_tga_file("./textures/fish.tga");
    texture5.flip_vertically(); // Flip the texture vertically because TGA files are stored upside down       
    
    texture6.read_tga_file("./textures/underwater.tga");
    texture6.flip_vertically(); // Flip the texture vertically because TGA files are stored upside down   

    //specular mapping for model
    specular.read_tga_file("./textures/specular/jellytexSpec.tga");
    specular.flip_vertically();

    specular2.read_tga_file("./textures/specular/seaweedSpec.tga");
    specular2.flip_vertically();

    specular3.read_tga_file("./textures/specular/sandSpec.tga");
    specular3.flip_vertically();

    specular4.read_tga_file("./textures/specular/coralSpec.tga");
    specular4.flip_vertically();

    specular5.read_tga_file("./textures/specular/fishSpec.tga");
    specular5.flip_vertically();

    specular6.read_tga_file("./textures/specular/underwaterSpec.tga");
    specular6.flip_vertically();
  
    // initialise SDL2
    init();

    //create vertexShader
    vertexShader vs;

    // Compute screen coordinates based on camera parameters
    float t, b, l, r;

    // Call a function to calculate screen coordinates based on the following parameters:
    // filmApertureWidth: width of the camera film aperture
    // filmApertureHeight: height of the camera film aperture
    // imageWidth: width of the image sensor or film plane
    // imageHeight: height of the image sensor or film plane
    // kOverscan: amount of overscan to apply to the image
    // nearClippingPlane: distance from the camera to the near clipping plane
    // focalLength: distance from the camera to the image plane
    // Store the results in t, b, l, and r
    computeScreenCoordinates(filmApertureWidth, filmApertureHeight, imageWidth, imageHeight, kOverscan, nearClippingPlane, focalLength, t, b, l, r);

    // Define the depth-buffer. Initialize depth buffer to far clipping plane.
    float* depthBuffer = new float[imageWidth * imageHeight]; 
    for (uint32_t i = 0; i < imageWidth * imageHeight; ++i) depthBuffer[i] = farClippingPlane;
    //define frame-buffer
    Vec3<unsigned char>* frameBuffer = new Vec3<unsigned char>[imageWidth * imageHeight];
    for (uint32_t i = 0; i < imageWidth * imageHeight; ++i) frameBuffer[i] = Vec3<unsigned char>(255);

    float camAngleX = 0.0f;
    float camAngleY = 0.0f;
    SDL_Event e;
    bool running = true;
    while (running) {
        // Start timer so we can gather frame generation statistics
        auto t_start = std::chrono::high_resolution_clock::now();

        // Clear back buffer, pixel data on surface and depth buffer (as movement), RGB sets the colour of the screen - in this case its dark turquoise 
        SDL_FillRect(screen, nullptr, SDL_MapRGB(screen->format, 4, 92, 90));
        SDL_RenderClear(renderer);
        // Only required if animating the camera as the depth buffer will need to be recomputed
        // Reset depth buffer to far clipping plane.
        for (uint32_t i = 0; i < imageWidth * imageHeight; ++i) depthBuffer[i] = farClippingPlane;


        //arcball camera
        float distance = 3.0f; //straight line distance between camera and look at point
        camAngleX++; //increment the x-axis camera angle
        if (camAngleX > 360.0f) camAngleX = 0.0f; //if camera angle is greater than 360, reset it to 0

        //calculate the camera's position via distance and angles
        float camX = distance * -sinf(camAngleX * (M_PI / 180)) * cosf((camAngleY) * (M_PI / 180));
        float camY = distance * -sinf((camAngleY) * (M_PI / 180));
        float camZ = distance * cosf((camAngleX) * (M_PI / 180)) * cosf((camAngleY) * (M_PI / 180));

        ////for cc_t model
        //Vec3f eye(camX, camY + 1, camZ); //create a vector representing the camera's position
        //Vec3f target(0.f, 0.0f, 0.f); //create a vector representing the camera's look-at point
        //Vec3f up(0.f, 1.f, 0.f); //create a vector representing the camera's up direction

        //for jelly
        Vec3f eye(-13.148364f, 26.358925f, 17.128357f);
        Vec3f target(5.305188f, 12.327177f, -14.828231f);
        Vec3f up(0.029893f, 0.997639f, -0.061837f);
        
        /*Vec3f eye(0.03f, 1.f, 30.f);
        Vec3f target(5.0f, 10.0f, 5.f);
        Vec3f up(0.f, 10.f, 0.f);*/

        worldToCamera = lookAt(eye, target, up).inverse(); //use the lookAt function to create a transformation matrix that maps world coordinates to camera coordinates, and then take its inverse to get the camera-to-world matrix

        Vec3f* vs0 = new Vec3f();
        Vec3f* vs1 = new Vec3f();
        Vec3f* vs2 = new Vec3f();

        TGAImage image(screen->w, screen->h, TGAImage::RGB);

        // create an array of models
        Model* models[6] = { model, model2, model3, model4, model5, model6 };
        TGAImage* textures[6] = { &texture, &texture2, &texture3, &texture4, &texture5, &texture6 };
        TGAImage* speculars[6] = { &specular, &specular2, &specular3, &specular4, &specular5, &specular6 };

        // iterate over each model in the array
        for (int m = 0; m < 6; m++) {
            Model* currentModel = models[m];
            TGAImage* currentTexture = textures[m];
            TGAImage* currentSpecular = speculars[m];

            for (uint32_t i = 0; i < currentModel->nfaces(); ++i) {

                // v0, v1 and v2 store the vertex positions of every vertex of the 3D model
                const Vec3f& v0 = currentModel->vert(currentModel->face(i)[0]);
                const Vec3f& v1 = currentModel->vert(currentModel->face(i)[1]);
                const Vec3f& v2 = currentModel->vert(currentModel->face(i)[2]);

                vs.processVertex(v0, vs0);
                vs.processVertex(v1, vs1);
                vs.processVertex(v2, vs2);

                // Convert the vertices of the triangle to raster space - you will need to implement convertToRaster()
                Vec3f v0Raster, v1Raster, v2Raster;
                convertToRaster(*vs0, worldToCamera, l, r, t, b, nearClippingPlane, imageWidth, imageHeight, v0Raster);
                convertToRaster(*vs1, worldToCamera, l, r, t, b, nearClippingPlane, imageWidth, imageHeight, v1Raster);
                convertToRaster(*vs2, worldToCamera, l, r, t, b, nearClippingPlane, imageWidth, imageHeight, v2Raster);

                // Precompute reciprocal of vertex z-coordinate
                v0Raster.z = 1 / v0Raster.z;  // reciprocal of v0 z-coordinate
                v1Raster.z = 1 / v1Raster.z;  // reciprocal of v1 z-coordinate
                v2Raster.z = 1 / v2Raster.z;  // reciprocal of v2 z-coordinate

                // Prepare vertex attributes. Divide them by their vertex z-coordinate
                // (though we use a multiplication here because v.z = 1 / v.z)
                // st0, st1 and st2 store the texture coordinates from the model of each vertex
                Vec2f st0 = currentModel->vt(currentModel->face(i)[0]);  // texture coordinates for v0
                Vec2f st1 = currentModel->vt(currentModel->face(i)[1]);  // texture coordinates for v1
                Vec2f st2 = currentModel->vt(currentModel->face(i)[2]);  // texture coordinates for v2

                Vec3f v0n = currentModel->vn(currentModel->vertex_normal(i)[0]);  //vertex normal for v0
                Vec3f v1n = currentModel->vn(currentModel->vertex_normal(i)[1]);  //vertex normal for v1
                Vec3f v2n = currentModel->vn(currentModel->vertex_normal(i)[2]);  //vertex normal for v2

                //this is needed for perspective correct interpolation
                st0 *= v0Raster.z;  // multiply st0 by reciprocal of v0 z-coordinate
                st1 *= v1Raster.z;  // multiply st1 by reciprocal of v1 z-coordinate
                st2 *= v2Raster.z;  // multiply st2 by reciprocal of v2 z-coordinate

                // Calculate the bounding box of the triangle defined by the vertices
                float xmin = min3(v0Raster.x, v1Raster.x, v2Raster.x);  // minimum x-coordinate of triangle
                float ymin = min3(v0Raster.y, v1Raster.y, v2Raster.y);  // minimum y-coordinate of triangle
                float xmax = max3(v0Raster.x, v1Raster.x, v2Raster.x);  // maximum x-coordinate of triangle
                float ymax = max3(v0Raster.y, v1Raster.y, v2Raster.y);  // maximum y-coordinate of triangle

                // the triangle is out of screen
                if (xmin > imageWidth - 1 || xmax < 0 || ymin > imageHeight - 1 || ymax < 0) continue;  // if the triangle is completely out of the screen, skip it and continue with the next one

                // sets the bounds of the rectangle for the raster triangle
                // be careful xmin/xmax/ymin/ymax can be negative. Don't cast to uint32_t
                uint32_t x0 = std::max(int32_t(0), (int32_t)(std::floor(xmin)));  // minimum x-coordinate of the bounding box
                uint32_t x1 = std::min(int32_t(imageWidth) - 1, (int32_t)(std::floor(xmax)));  // maximum x-coordinate of the bounding box
                uint32_t y0 = std::max(int32_t(0), (int32_t)(std::floor(ymin)));  // minimum y-coordinate of the bounding box
                uint32_t y1 = std::min(int32_t(imageHeight) - 1, (int32_t)(std::floor(ymax)));

                // calculates the area of the triangle, used in determining barycentric coordinates
                float area = edgeFunction(v0Raster, v1Raster, v2Raster);

                for (uint32_t y = y0; y <= y1; ++y) {
                    for (uint32_t x = x0; x <= x1; ++x) {

                        // Initialize accumulators for color and depth values
                        Vec3f accumColour(0, 0, 0);
                        float accumDepth = 0;
                        int sampleCount = 0;

                        // Sample 4x4 grid for supersampling anti aliasing
                        for (float dy = 0; dy < 1; dy += 0.25) {
                            for (float dx = 0; dx < 1; dx += 0.25) {
                                // Sample a pixel at the current position within the supersampling grid
                                Vec3f pixelSample(x + dx, y + dy, 0);

                                // Calculate the area of the subtriangles for barycentric coordinates
                                float w0 = edgeFunction(v1Raster, v2Raster, pixelSample);
                                float w1 = edgeFunction(v2Raster, v0Raster, pixelSample);
                                float w2 = edgeFunction(v0Raster, v1Raster, pixelSample);

                                if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
                                    // divide by the area to give us our coefficients
                                    w0 /= area;
                                    w1 /= area;
                                    w2 /= area;

                                    // Interpolate z values using barycentric coordinates
                                    float oneOverZ = v0Raster.z * w0 + v1Raster.z * w1 + v2Raster.z * w2;
                                    float z = 1 / oneOverZ;

                                    // Depth-buffer test
                                    if (z < depthBuffer[y * imageWidth + x]) { // Is this triangle closer than others previously? 
                                        //this makes sure that only triangles closer to the camera get drawn
                                        depthBuffer[y * imageWidth + x] = z;  // Update the depth buffer

                                        // Calculate the texture coordinate based on barycentric position of the pixel
                                        Vec2f st = st0 * w0 + st1 * w1 + st2 * w2;

                                        // correct for perspective distortion
                                        st *= z;

                                        //assigns textureColour so not hardcoded to one texture colour
                                        TGAColor textureColour = currentTexture->get(st.x * currentTexture->get_width(), st.y * currentTexture->get_height()); //sample texture colour at computed coord st
                                        //assigns specularColour
                                        TGAColor specularColour = currentSpecular->get(st.x * currentSpecular->get_width(), st.y * currentSpecular->get_height());

                                        // If need to compute the actual position of the shaded
                                        // point in camera space. Proceed like with the other vertex attribute.
                                        // Divide the point coordinates by the vertex z-coordinate then
                                        // interpolate using barycentric coordinates and finally multiply
                                        // by sample depth.
                                        Vec3f v0Cam, v1Cam, v2Cam;
                                        worldToCamera.multVecMatrix(v0, v0Cam);
                                        worldToCamera.multVecMatrix(v1, v1Cam);
                                        worldToCamera.multVecMatrix(v2, v2Cam);

                                        //divide them by the respective z-coordinate as with any other vertex attribute and interpolate using barycentric coordinates
                                        float px = (v0Cam.x / -v0Cam.z) * w0 + (v1Cam.x / -v1Cam.z) * w1 + (v2Cam.x / -v2Cam.z) * w2;
                                        float py = (v0Cam.y / -v0Cam.z) * w0 + (v1Cam.y / -v1Cam.z) * w1 + (v2Cam.y / -v2Cam.z) * w2;

                                        //p in camera space
                                        Vec3f pt(px * z, py * z, -z); // pt is in camera space

                                        // Compute the face normal which is used for a simple facing ratio.
                                        // Keep in mind that we are doing all calculation in camera space.
                                        // Thus the view direction can be computed as the point on the object
                                        // in camera space minus Vec3f(0), the position of the camera in camera space.
                                        Vec3f n = (v1Cam - v0Cam).crossProduct(v2Cam - v0Cam);
                                        n.normalize();
                                        Vec3f viewDirection = -pt;
                                        viewDirection.normalize();

                                        // Calculate shading of the surface based on dot product of the normal and view direction
                                        float nDotView = std::max(0.f, n.dotProduct(viewDirection));

                                        // Calculate the shading for each vertex and interpolate to find shading for the current pixel
                                        // Interpolation coefficients w0, w1, and w2 are precomputed
                                        Vec3f sn = v0n * w0 + v1n * w1 + v2n * w2;
                                        nDotView = std::max(0.f, sn.dotProduct(viewDirection));

                                        // Now to calculate specular reflection
                                        // Compute the half-vector
                                        Vec3f textureColorNormalized = Vec3f(textureColour.r / 255.0f, textureColour.g / 255.0f, textureColour.b / 255.0f);
                                        Vec3f diffuseComponent = nDotView * textureColorNormalized;

                                        // Now to calculate specular reflection
                                        // Compute the half-vector
                                        Vec3f halfVector = (viewDirection + lightDirection).normalize();

                                        // Calculate the dot product of the half-vector and the normal
                                        float nDotH = std::max(0.f, sn.dotProduct(halfVector));

                                        // Calculate the specular component
                                        float specular = pow(nDotH, shininess);
                                        Vec3f specularComponent = specular * lightColour;

                                        // Display only the diffuse component
                                        //Vec3f finalColor = diffuseComponent;

                                        // Uncomment the following line to display attempted the specular mapping
                                        Vec3f finalColor = specularComponent + diffuseComponent;

                                        float red = finalColor.x;
                                        float green = finalColor.y;
                                        float blue = finalColor.z;

                                        // The final color is the result of the fraction multiplied by the
                                        // checkerboard pattern defined in checker.
                                       /* const int M = 10;
                                        float checker = (fmod(st.x * M, 1.0) > 0.5) ^ (fmod(st.y * M, 1.0) < 0.5);
                                        float c = 0.3 * (1 - checker) + 0.7 * checker;*/
                                        //nDotView *= c;

                                        //// Set the x component of the pixel at coordinates (x, y) in the frameBuffer to red times 255.
                                        //frameBuffer[y * imageWidth + x].x = red * 255;

                                        //// Set the y component of the pixel at coordinates (x, y) in the frameBuffer to green times 255.
                                        //frameBuffer[y * imageWidth + x].y = green * 255;

                                        //// Set the z component of the pixel at coordinates (x, y) in the frameBuffer to blue times 255.
                                        //frameBuffer[y * imageWidth + x].z = blue * 255;

                                        //image.set(x, y, TGAColor((unsigned char)(red * 255), (unsigned char)(green * 255), (unsigned char)(blue * 255), 255));

                                        //// Set the pixel value on the SDL_Surface that gets drawn to the SDL_Window
                                        //Uint32 colour = SDL_MapRGB(screen->format, red * 255, green * 255, blue * 255);
                                        //putpixel(screen, x, y, colour);

                                         // Accumulate color and depth values
                                        accumColour.x += red * 255;
                                        accumColour.y += green * 255;
                                        accumColour.z += blue * 255;
                                        accumDepth += z;
                                        sampleCount++;
                                    }
                                }
                                // If any samples passed the depth test, calculate the average color and depth and set the pixel value
                                if (sampleCount > 0) {
                                    // Calculate average color and depth values
                                    Vec3f avgColour(accumColour.x / sampleCount, accumColour.y / sampleCount, accumColour.z / sampleCount);
                                    float avgDepth = accumDepth / sampleCount;

                                    // Update the depth buffer with the average depth
                                    depthBuffer[y * imageWidth + x] = avgDepth;

                                    // Set the pixel value in the frameBuffer and image
                                    frameBuffer[y * imageWidth + x].x = avgColour.x;
                                    frameBuffer[y * imageWidth + x].y = avgColour.y;
                                    frameBuffer[y * imageWidth + x].z = avgColour.z;
                                    image.set(x, y, TGAColor((unsigned char)(avgColour.x), (unsigned char)(avgColour.y), (unsigned char)(avgColour.z), 255));

                                    // Set the pixel value on the SDL_Surface that gets drawn to the SDL_Window
                                    Uint32 colour = SDL_MapRGB(screen->format, avgColour.x, avgColour.y, avgColour.z);
                                    putpixel(screen, x, y, colour);
                                }
                            }
                        }
                    }
                }
            }
        }
        vs0 = nullptr; vs1 = nullptr; vs2 = nullptr;
        delete vs0, vs1, vs2;

        image.write_tga_file("rasteriser.tga");

        // Calculate frame interval timing
        auto t_end = std::chrono::high_resolution_clock::now();
        auto passedTime = std::chrono::duration<double, std::milli>(t_end - t_start).count();
        std::cerr << "Frame render time: " << passedTime << " ms" << std::endl;
        
        // Create texture from the surface and RenderCopy/Present from backbuffer
        SDL_Texture* texture = SDL_CreateTextureFromSurface(renderer, screen);
        if (texture == NULL) {
            fprintf(stderr, "CreateTextureFromSurface failed: %s\n", SDL_GetError());
            exit(1);
        }
        SDL_FreeSurface(screen);

        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);

        // Clean up heap allocation
        SDL_DestroyTexture(texture);

        // Check for ESC sequence, otherwise keep drawing frames
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

    // tidy up dangling pointer to the depth-buffer
    delete[] depthBuffer;
    //tidy up frame-buffer
    delete[] frameBuffer;

    delete model, model2, model3, model4, model5;

    return 0;
}