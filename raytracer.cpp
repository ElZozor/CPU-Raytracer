#include "image.h"
#include "kdtree.h"
#include "ray.h"
#include "raytracer.h"
#include "scene_types.h"
#include <stdio.h>
#include <cmath>
#include <map>
#include <algorithm>
#include <cstring>
#include <functional>


#include "intersections.h"

#include <glm/gtc/epsilon.hpp>

#define BETWEEN3D(position, minimum, maximum)    (position.x >= minimum.x && position.x <= maximum.x        \
                                                    && position.y >= minimum.y && position.y <= maximum.y   \
                                                    && position.z >= minimum.z && position.z <= maximum.z)



/// acne_eps is a small constant used to prevent acne when computing
/// intersection
//    or boucing (add this amount to the position before casting a new ray !
const float acne_eps = 1e-4;



// Custom variables that are used to
// perform the depth of field filter
static std::vector<float> heightMap;
std::string SceneParameters::AAMode     = "sobel";
float SceneParameters::focalDistance    = -1.f;
float SceneParameters::focalRange       = 2.f;
int SceneParameters::antiAliasing       = 4;
size_t SceneParameters::imageWidth      = 800;
size_t SceneParameters::imageHeight     = 600;


float gaussianCoeff(size_t x, size_t y, int matrixSize);

/* ---------------------------------------------------------------------------
 */
/*
 *	The following functions are coded from Cook-Torrance bsdf model
 *description and are suitable only
 *    for rough dielectrics material (RDM. Code has been validated with Mitsuba
 *renderer)
 */

// Shadowing and masking function. Linked with the NDF. Here, Smith function,
// suitable for Beckmann NDF
float RDM_chiplus(float c) { return (c > 0.f) ? 1.f : 0.f; }

/** Normal Distribution Function : Beckmann
 * NdotH : Norm . Half
 */

float RDM_Beckmann(float NdotH, float alpha)
{
    const float cosSquared   = NdotH * NdotH;
    const float tanOHSquared = (1 - cosSquared) / cosSquared;
    const float alphaSquared = alpha * alpha;
    const float numerateur   = expf((-tanOHSquared) / (alphaSquared));
    const float denominateur = M_PI * alphaSquared * (cosSquared * cosSquared);

    
    return RDM_chiplus(NdotH) * (numerateur / denominateur);
}

// Fresnel term computation. Implantation of the exact computation. we can use
// the Schlick approximation
// LdotH : Light . Half
float RDM_Fresnel(float cosOi, float n1, float n2)
{
    const float n1dn2 = (n1 / n2);
    float sin2Ot = (n1dn2 * n1dn2) * (1 - (cosOi * cosOi));
    if (sin2Ot > 1.f)
    {
        return 1.f;
    }

    const float cosOt = sqrt(1.f - sin2Ot);

    const float rs = (powf(n1 * cosOi - n2 * cosOt, 2)) / (powf(n1 * cosOi + n2 * cosOt, 2));
    const float rp = (powf(n1 * cosOt - n2 * cosOi, 2)) / (powf(n1 * cosOt + n2 * cosOi, 2));

    return 0.5f * (rs + rp);
}

// DdotH : Dir . Half
// HdotN : Half . Norm
float RDM_G1(float DdotH, float DdotN, float alpha)
{
    const float tanOx = (sqrtf(1.f - (DdotN * DdotN)) / DdotN);
    const float b = (1 / (alpha * tanOx));
    const float k = (DdotH / DdotN);

    if (b < 1.6f)
    {
        return RDM_chiplus(k) * ((3.535f * b + 2.181f * (b * b)) / (1.f + 2.276f * b + 2.577f * (b * b)));
    }
    else
    {
        return RDM_chiplus(k);
    }
}

// LdotH : Light . Half | v
// LdotN : Light . Norm | l
// VdotH : View . Half    | h
// VdotN : View . Norm    | n
float RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN,
                                float alpha)
{

    const float G1A = RDM_G1(LdotH, LdotN, alpha);
    const float G1B = RDM_G1(VdotH, VdotN, alpha);

    return G1A * G1B;
}

// Specular term of the Cook-torrance bsdf
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN,
                                    float VdotN, Material *m)
{
    const auto& ks =  m->specularColor;
    const float D = RDM_Beckmann(NdotH, m->roughness);
    const auto F = RDM_Fresnel(LdotH, 1.f, m->IOR);
    const auto G = RDM_Smith(LdotH, LdotN, VdotH, VdotN, m->roughness);

    return ks * ((D*F*G) / (4.f * LdotN * VdotN));
}

// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(Material *m)
{
    return m->diffuseColor / M_PIf32;
}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN,
                                float VdotN, Material *m)
{
    const auto rightTerm = RDM_bsdf_s(LdotH, NdotH, VdotH, LdotN, VdotN, m);
    const auto leftTerm = RDM_bsdf_d(m);

    return (leftTerm + rightTerm);
}

color3 shade(vec3 n, vec3 v, vec3 l, color3 lc, Material *mat)
{
    const auto h = normalize(v+l);
    const auto LdotH = dot(l, h), NdotH = dot(n, h),
        VdotH = dot(v, h), LdotN = dot(l, n), VdotN = dot(v, n);
    
    color3 ret = lc * RDM_bsdf(LdotH, NdotH, VdotH, LdotN, VdotN, mat) * LdotN;
    return clamp(ret, vec3(0.f), vec3(1.f));
}

//! if tree is not null, use intersectKdTree to compute the intersection instead
//! of intersect scene


color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree)
{
    if (ray->depth >= 10)
    {
        return color3(0);
    }

    Intersection intersection;
    if (!intersectKdTree(scene, tree, ray, &intersection))
    {
        return scene->skyColor;
    }


    color3 color(0);
    for (const auto &light : scene->lights)
    {
        const auto &L = light->position;
        const auto &P = intersection.position;
        const auto l = normalize(L - P);

        Ray shadowRay;
        rayInit(&shadowRay, P, l, acne_eps, distance(L, P));
        Intersection shadowIntersection;

        if (!intersectObjectKdTree(scene, tree, &shadowRay, &shadowIntersection))
        {
            color += shade(intersection.normal, -ray->dir, l, light->color, intersection.mat);
        }
    }


    const auto reflectionDir = reflect(ray->dir, intersection.normal);
    const int add = max(1, 10 - int(intersection.mat->IOR * 10));
    
    Ray reflectionRay;
    rayInit(&reflectionRay, intersection.position + (reflectionDir * 0.001f), reflectionDir, acne_eps, 100000.f, ray->depth + add);

    const color3 cr = trace_ray(scene, &reflectionRay, tree);
    const float F = min(1.f, RDM_Fresnel(dot(reflectionRay.dir, intersection.normal), 1.f, intersection.mat->IOR));

    return color + F * cr * intersection.mat->specularColor;
}


// /*
//    _____                           _                       _             
//   / ____|                         | |                     (_)            
//  | (___  _   _ _ __   ___ _ __ ___| | __ _ _ __ ___  _ __  _ _ __   __ _ 
//   \___ \| | | | '_ \ / _ \ '__/ __| |/ _` | '_ ` _ \| '_ \| | '_ \ / _` |
//   ____) | |_| | |_) |  __/ |  \__ \ | (_| | | | | | | |_) | | | | | (_| |
//  |_____/ \__,_| .__/ \___|_|  |___/_|\__,_|_| |_| |_| .__/|_|_| |_|\__, |
//               | |                                   | |             __/ |
//               |_|                                   |_|            |___/ 

// */

/**
 * @brief Naïve version of AA.
 * 
 * @param scene         The scene to render
 * @param tree          The kdtree to compute code
 * @param ray_delta_x   
 * @param ray_delta_y 
 * @param x 
 * @param y 
 * @param dx 
 * @param dy 
 * @param amount 
 * @return color3 
 */
color3 superslamping(Scene* scene, KdTree* tree, const vec3& ray_delta_x, const vec3& ray_delta_y, float x, 
                    float y, const vec3& dx, const vec3& dy) {
    
    const auto AA = SceneParameters::antiAliasing;
    color3 result(0);
    for (int xx = 0; xx < AA; ++xx)
    {
        for (int yy = 0; yy < AA; ++yy)
        {
            vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
                (x + (static_cast<float>(xx) / AA)) * dx + (y + (static_cast<float>(yy) / AA)) * dy;

            Ray rx;
            rayInit(&rx, scene->cam.position, normalize(ray_dir));

            result += trace_ray(scene, &rx, tree);
        }
    }

    return result / float(AA * AA);

}


/*
   _____       _         _ _       _     _             
  / ____|     | |       | (_)     (_)   (_)            
 | (___  _   _| |__   __| |___   ___ ___ _  ___  _ __  
  \___ \| | | | '_ \ / _` | \ \ / / / __| |/ _ \| '_ \ 
  ____) | |_| | |_) | (_| | |\ V /| \__ \ | (_) | | | |
 |_____/ \__,_|_.__/ \__,_|_| \_/ |_|___/_|\___/|_| |_|
                                                       
*/                                                


/**
 * @brief Subdivide a pixel, based on the AA amount specified.
 * 
 * @param img           The resulting image
 * @param y             The y position of the pixel
 * @param x             The x position of the pixel
 * @param scene         The scene which contains the objects
 * @param tree          The KDTree used to compute colisions
 * @param ray_delta_x   The dx of the ray
 * @param ray_delta_y   The dy of the ray
 * @param dx            The dx
 * @param dy            The dy
 * @return color3       The resulting pixel
 */
color3 subdividePixel(Image* img, size_t y, size_t x, Scene* scene, KdTree* tree, 
                        const vec3& ray_delta_x, const vec3& ray_delta_y,
                        const vec3& dx, const vec3& dy)
{
    // First we test if the AA paramater is
    // greater than 1, otherwise we simply return
    // the already computed pixel.
    const auto AA = SceneParameters::antiAliasing;
    if (AA < 2)
    {
        return img->data[y * img->width + x];
    }


    // Then we subdivide the pixel into subpixels, 
    // compute the resulting value for each
    // of them and return the average value.
    color3 result(img->data[y * img->width + x]);
    for (int xx = 0; xx < AA; ++xx)
    {
        for (int yy = 0; yy < AA; ++yy)
        {
            if (xx + yy == 0) continue;
            
            Ray rx;
            vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
                (x + (static_cast<float>(xx) / AA)) * dx + (y + (static_cast<float>(yy) / AA)) * dy;

            rayInit(&rx, scene->cam.position, normalize(ray_dir));

            result += trace_ray(scene, &rx, tree);
        }
    }

    return result / (float(AA * AA));
}



/*
  _                           
 | |                          
 | |    _   _ _ __ ___   __ _ 
 | |   | | | | '_ ` _ \ / _` |
 | |___| |_| | | | | | | (_| |
 |______\__,_|_| |_| |_|\__,_|
                              
*/

// Compute the luma value
#define rgb2luma(rgbcolor) (sqrtf(dot(rgbcolor, vec3(0.299, 0.587, 0.114))))

// Construct the luma array
#define CONSTRUCT_ARRAY(img, y, x) {    \
    rgb2luma(img->data[(y-1) * img->width + (x-1)]),                    \
    rgb2luma(img->data[(y-1) * img->width + (x+0)]),                    \
    rgb2luma(img->data[(y-1) * img->width + (x+1)]),                    \
    rgb2luma(img->data[(y+0) * img->width + (x-1)]),                    \
    rgb2luma(img->data[(y+0) * img->width + (x+0)]),                    \
    rgb2luma(img->data[(y+0) * img->width + (x+1)]),                    \
    rgb2luma(img->data[(y+1) * img->width + (x-1)]),                    \
    rgb2luma(img->data[(y+1) * img->width + (x+0)]),                    \
    rgb2luma(img->data[(y+1) * img->width + (x+1)])                     \
}

/**
 * @brief Improve render time when applying AA to the \
 * final render. The AA is only applied to the edges.
 * 
 * @param img           The resulting image
 * @param y             The y position of the pixel
 * @param x             The x position of the pixel
 * @param scene         The scene which contains the objects
 * @param tree          The KDTree used to compute colisions
 * @param ray_delta_x   The dx of the ray
 * @param ray_delta_y   The dy of the ray
 * @param dx            The dx
 * @param dy            The dy
 * @return color3       The resulting pixel
 */
color3 computePixelLumaMethod(Image* img, size_t y, size_t x, Scene* scene, KdTree* tree, 
                                    const vec3& ray_delta_x, const vec3& ray_delta_y,
                                    const vec3& dx, const vec3& dy)
{
    // If we are on the edge of the picture, we simply
    // return the pixel color.
    const auto& pixelColor = img->data[y * img->width + x];
    if (y == 0 || y == (img->height - 1))   return pixelColor;
    if (x == 0 || x == (img->width - 1))    return pixelColor;

    // Here we compute the luma array for the given pixel
    // and it's neighbors. We then find the minimum and
    // maximum value to compute the range between them.
    float lumaCoeff[9] = CONSTRUCT_ARRAY(img, y, x);
    float minLuma = lumaCoeff[0], maxLuma = lumaCoeff[0];
    for (size_t i = 1; i < 9; ++i)
    {
        minLuma = std::min(minLuma, lumaCoeff[i]);
        maxLuma = std::max(maxLuma, lumaCoeff[i]);
    }

    // Finally we test if we are under a certain 
    // thresold, in this case we return the pixel color
    // otherwise, we apply the AA.
    const float lumaRange = maxLuma - minLuma + 0.05f;
    if (lumaRange < std::max(0.0312f, maxLuma * 0.125f))
    {
        return pixelColor;
    }

    return subdividePixel(img, y, x, scene, tree, ray_delta_x, ray_delta_y, dx, dy);
}



/**
 * @brief This function simply iterate through the pixels \
 * call the luma function to compute pixels values.
 * 
 * @param img            The resulting image
 * @param scene         The scene which contains the objects
 * @param tree          The KDTree used to compute colisions
 * @param ray_delta_x   The dx of the ray
 * @param ray_delta_y   The dy of the ray
 * @param dx            The dx
 * @param dy            The dy
 */
void lumaImprovedSuperslamping(Image* img, Scene* scene, KdTree* tree, const vec3& ray_delta_x, const vec3& ray_delta_y,
                        const vec3& dx, const vec3& dy, int amount) {
    
    const size_t height = SceneParameters::imageHeight;
    const size_t width  = SceneParameters::imageWidth;
    std::vector<color3> buffer(height*width);

    #pragma omp parallel for schedule(dynamic)
    for (size_t y = 0; y < height; ++y)
    {
        for (size_t x = 0; x < width; ++x)
        {
            buffer[y * img->width + x] = computePixelLumaMethod(img, y, x, scene, tree, ray_delta_x, ray_delta_y, dx, dy);
        }
    }

    memcpy(img->data, buffer.data(), sizeof(color3) * height * width);

}



/*
   _____       _          _ 
  / ____|     | |        | |
 | (___   ___ | |__   ___| |
  \___ \ / _ \| '_ \ / _ \ |
  ____) | (_) | |_) |  __/ |
 |_____/ \___/|_.__/ \___|_|

*/

/**
 * @brief Detect the horizontal edges for  \
 * the sobel operator
 * 
 */
#define SOBEL_X_VALUE(buffer, y, x) (               \
      buffer[(y - 1) * img->width + (x - 1)]        \
    + buffer[(y + 0) * img->width + (x - 1)] * (+2) \
    + buffer[(y + 1) * img->width + (x - 1)]        \
    + buffer[(y - 1) * img->width + (x + 1)] * (-1) \
    + buffer[(y + 0) * img->width + (x + 1)] * (-2) \
    + buffer[(y + 1) * img->width + (x + 1)] * (-1) \
)

/**
 * @brief Detect the vertical edges for  \
 * the sobel operator
 * 
 */
#define SOBEL_Y_VALUE(buffer, y, x) (               \
      buffer[(y - 1) * img->width + (x - 1)]        \
    + buffer[(y - 1) * img->width + (x + 0)] * (+2) \
    + buffer[(y - 1) * img->width + (x + 1)]        \
    + buffer[(y + 1) * img->width + (x - 1)] * (-1) \
    + buffer[(y + 1) * img->width + (x + 0)] * (-2) \
    + buffer[(y + 1) * img->width + (x + 1)] * (-1) \
)


/**
 * @brief Mix both horizontal and vertical values
 * 
 */
#define SOBEL_VALUE(img, y, x) (                    \
      abs(SOBEL_X_VALUE(img, y, x))                 \
    + abs(SOBEL_Y_VALUE(img, y, x))                 \
)

/**
 * @brief Transform a pixel into it's greyscale value
 * 
 */
#define GREY_SCALE(color) (                         \
      color[0] * 0.07f                              \
    + color[1] * 0.72f                              \
    + color[2] * 0.21f                              \
)



/**
 * @brief Improve the render time by detecting edges for superslamping
 * 
 * @param img           The resulting image
 * @param scene         The scene which contains the objects
 * @param tree          The KDTree used to compute colisions
 * @param ray_delta_x   The dx of the ray
 * @param ray_delta_y   The dy of the ray
 * @param dx            The dx
 * @param dy            The dy
 */
void sobelImprovedSuperslamping(Image* img, Scene* scene, KdTree* tree, 
                        const vec3& ray_delta_x, const vec3& ray_delta_y,
                        const vec3& dx, const vec3& dy)
{
    const size_t width = img->width, height = img->height;
    std::vector<int> greyscaleImage(width * height);
    std::vector<color3> sobelImage(width * height);

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < width * height; ++i)
    {
        greyscaleImage[i] = (GREY_SCALE(img->data[i]) * 255);
    }

    for (size_t i = width; i < height * width - width; i+=width)
    {
        sobelImage[i] = img->data[i];
        sobelImage[i + (width - 1)] = img->data[i + (width - 1)];
    }
    std::memcpy(sobelImage.data(), img->data, sizeof(color3) * width);
    std::memcpy(sobelImage.data() + (height - 1) * width, img->data + (height - 1) * width, sizeof(color3) * width);

    #pragma omp parallel for schedule(dynamic)
    for (size_t y = 1; y < height - 1; ++y)
    {
        for (size_t x = 1; x < width - 1; ++x)
        {
            if (SOBEL_VALUE(greyscaleImage, y, x) > 0)
            {
                sobelImage[y * width + x] = subdividePixel(img, y, x, scene, tree, ray_delta_x, ray_delta_y, dx, dy);
            }
            else
            {
                sobelImage[y * width + x] = img->data[y * width + x];
            }
        }
    }

    std::memcpy(img->data, sobelImage.data(), sizeof(color3) * width * height);
}



/**
 * @brief Compute the heightmap for the gaussian blur
 * 
 * @param scene         The scene which contains the objects
 * @param tree          The KDTree used to compute colisions
 * @param ray_delta_x   The dx of the ray
 * @param ray_delta_y   The dy of the ray
 * @param dx            The dx
 * @param dy            The dy
 */
void computeHeightMap(Scene* scene, KdTree* tree, const vec3& ray_delta_x, 
        const vec3& ray_delta_y,const vec3& dx, const vec3& dy) 
{

    const size_t width  = SceneParameters::imageWidth;
    const size_t height = SceneParameters::imageHeight;
    heightMap.reserve(width * height);

    // The schedule is set to dynamic due to the different
    // render times for different pixels.
    #pragma omp parallel for schedule(dynamic)
    for (size_t y = 0; y < height; ++y)
    {
        for (size_t x = 0; x < width; ++x)
        {
            vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
                float(x) * dx + float(y) * dy;

            Ray rx;
            rayInit(&rx, scene->cam.position, normalize(ray_dir));
            Intersection intersection;

            if (intersectKdTree(scene, tree, &rx, &intersection))
            {
                const float distance = rx.tmax;
                heightMap[y * width + x] = min(8.f, max(0.f, abs(SceneParameters::focalDistance - distance) - SceneParameters::focalRange));
            }
            else
            {
                heightMap[y * width + x] = 8.f;
            }
        }
    }


}




/*
  ____  _            
 |  _ \| |           
 | |_) | |_   _ _ __ 
 |  _ <| | | | | '__|
 | |_) | | |_| | |   
 |____/|_|\__,_|_|   
                     
*/


/**
 * @brief Compute the pascal triangle given line
 * 
 * @param size                  The desired line
 * @return std::vector<int>     The resulting line
 */
std::vector<int> pascal(int size)
{
    if (size <= 1)
    {
        return std::vector<int>({1});
    }

    const auto p = pascal(size - 1);
    const auto psize = p.size();
    std::vector<int> result({p[0]});
    for (int i = 0; i < psize - 1; ++i)
    {
        result.push_back(p[i] + p[i + 1]);
    }
    result.push_back(p[psize - 1]);

    return result;
}



/**
 * @brief Approximate the gaussian coeff \
 *  via the pascal triangle
 * 
 * @param x             The x position of the pixel
 * @param y             The y position of the pixel
 * @param matrixSize    The desired matrix size
 * @return float        The approximed coefficient
 */
float gaussianCoeff(size_t x, size_t y, int matrixSize)
{
    // Here the variables are declared static to the function
    // to improve the render time.
    // We compute them one time ( has the maximum value is set to 16 )
    // and then each time we simply return the asked value.
    static bool initialized = false;
    static std::vector<int> pascalTriangle[16];

    if (initialized == false)
    {
        for (int i = 0; i < 16; ++i)
        {
            pascalTriangle[i] = pascal(i + 1);
        }

        initialized = true;
    }

    matrixSize = max(1, matrixSize);
    const auto& p = pascalTriangle[matrixSize - 1];
    
    return float(p[x] * p[y]);
}



/**
 * @brief Apply a gaussian blur to the image based \
 * on a computed heightmap.
 * 
 * @param img           The resulting image
 * @param scene         The scene which contains the objects
 * @param tree          The KDTree used to compute colisions
 * @param ray_delta_x   The dx of the ray
 * @param ray_delta_y   The dy of the ray
 * @param dx            The dx
 * @param dy            The dy
 */
void gaussianblur(Image* img, Scene* scene, KdTree* tree, const vec3& ray_delta_x, 
        const vec3& ray_delta_y,const vec3& dx, const vec3& dy)
{
    computeHeightMap(scene, tree, ray_delta_x, ray_delta_y, dx, dy);
    std::vector<color3> blurredImage(img->width * img->height);

    for (size_t y = 0; y < img->height; ++y)
    {
        for (size_t x = 0; x < img->width; ++x)
        {
            long offset = long(heightMap[y * img->width + x]);


            const size_t yymin = max(0L, long(y) - offset);
            const size_t yymax = min(img->height, y + offset);
            const size_t xxmin = max(0L, long(x) - offset);
            const size_t xxmax = min(img->width, x + offset);
            const size_t matrixSize = offset * 2;


            color3 accumulator(0);
            float compt(0);
            for (size_t yy = yymin; yy < yymax; ++yy)
            {
                const auto position = yy * img->width;
                for (size_t xx = xxmin; xx < xxmax; ++xx)
                {
                    if (heightMap[position + xx] > 0.f)
                    {
                        const float gaussian_coeff = gaussianCoeff((xx - xxmin), (yy - yymin), matrixSize);

                        accumulator += img->data[position + xx] * gaussian_coeff;
                        compt += gaussian_coeff;
                    }
                }
            }

            if (compt == 0.f)
            {
                blurredImage[y * img->width + x] = img->data[y * img->width + x];
            }
            else
            {
                blurredImage[y * img->width + x] = accumulator / compt;
            }
        }
    }

    memcpy(img->data, blurredImage.data(), sizeof(color3) * img->height * img->width);
}





void classicRender(Image* img, Scene* scene)
{
    float aspect = 1.f / scene->cam.aspect;
    KdTree *tree = initKdTree(scene);

    float delta_y = 1.f / (img->height * 0.5f);     //! one pixel size
    vec3 dy = delta_y * aspect * scene->cam.ydir; //! one pixel step
    vec3 ray_delta_y = (0.5f - img->height * 0.5f) / (img->height * 0.5f) *
                                         aspect * scene->cam.ydir;

    float delta_x = 1.f / (img->width * 0.5f);
    vec3 dx = delta_x * scene->cam.xdir;
    vec3 ray_delta_x =
            (0.5f - img->width * 0.5f) / (img->width * 0.5f) * scene->cam.xdir;

    size_t percent = 0;
    const size_t percent_range = img->height / 100.f;

    for (size_t j = 0; j < img->height; j++)
    {
        if (j > percent)
        {   
            if (j != 0)
                printf("\033[A\r");
            float progress = (float)j / img->height * 100.f;
            printf("progress\t[");
            int cpt = 0;
            for (cpt = 0; cpt < progress; cpt += 5)
                printf(".");
            for (; cpt < 100; cpt += 5)
                printf(" ");
            printf("]\n");
            percent = j + percent_range;
        }
        
#pragma omp parallel for
        for (size_t i = 0; i < img->width; i++)
        {
            color3 *ptr = getPixelPtr(img, i, j);

            vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
                                            float(i) * dx + float(j) * dy;

            Ray rx;
            rayInit(&rx, scene->cam.position, normalize(ray_dir));
            *ptr = trace_ray(scene, &rx, tree);
        }
    }


    if (SceneParameters::antiAliasing > 1)
    {
        printf("Computing Anti aliasing..\n");
        if (SceneParameters::AAMode == "luma")
        {
            lumaImprovedSuperslamping(img, scene, tree, ray_delta_x, ray_delta_y,
                                        dx, dy, SceneParameters::antiAliasing);
        }
        else
        {
            sobelImprovedSuperslamping(img, scene, tree, ray_delta_x, ray_delta_y,
                                        dx, dy);
        }
    }
    


    if (SceneParameters::focalDistance > 0)
    {
        printf("Computing depth of field..\n");
        gaussianblur(img, scene, tree, ray_delta_x, ray_delta_y, dx, dy);
    }
}



void naiveSuperslampingRender(Image* img, Scene* scene)
{
    float aspect = 1.f / scene->cam.aspect;
    KdTree *tree = initKdTree(scene);

    float delta_y = 1.f / (img->height * 0.5f);     //! one pixel size
    vec3 dy = delta_y * aspect * scene->cam.ydir; //! one pixel step
    vec3 ray_delta_y = (0.5f - img->height * 0.5f) / (img->height * 0.5f) *
                                         aspect * scene->cam.ydir;

    float delta_x = 1.f / (img->width * 0.5f);
    vec3 dx = delta_x * scene->cam.xdir;
    vec3 ray_delta_x =
            (0.5f - img->width * 0.5f) / (img->width * 0.5f) * scene->cam.xdir;

    size_t percent = 0;
    const size_t percent_range = img->height / 100.f;

    for (size_t j = 0; j < img->height; j++)
    {
        if (j > percent)
        {   
            if (j != 0)
                printf("\033[A\r");
            float progress = (float)j / img->height * 100.f;
            printf("progress\t[");
            int cpt = 0;
            for (cpt = 0; cpt < progress; cpt += 5)
                printf(".");
            for (; cpt < 100; cpt += 5)
                printf(" ");
            printf("]\n");
            percent = j + percent_range;
        }
        
#pragma omp parallel for
        for (size_t i = 0; i < img->width; i++)
        {
            color3 *ptr = getPixelPtr(img, i, j);
            *ptr = superslamping(scene, tree, ray_delta_x, ray_delta_y, i, j, dx, dy);
        }
    }

    if (SceneParameters::focalDistance > 0)
    {
        printf("Computing depth of field..\n");
        gaussianblur(img, scene, tree, ray_delta_x, ray_delta_y, dx, dy);
    }
}





void renderImage(Image *img, Scene *scene)
{    
    if (SceneParameters::AAMode == "none")
    {
        naiveSuperslampingRender(img, scene);
    }
    else
    {
        classicRender(img, scene);
    }
}
