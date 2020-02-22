#include "image.h"
#include "kdtree.h"
#include "ray.h"
#include "raytracer.h"
#include "scene_types.h"
#include <stdio.h>
#include <cmath>
#include <map>


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
float SceneParameters::focalDistance = -1.f;
float SceneParameters::focalRange = 2.f;
int SceneParameters::antiAliasing = 4;
size_t SceneParameters::imageWidth;
size_t SceneParameters::imageHeight;



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
    return m->diffuseColor / (float)M_PI;
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

    // const auto bsdf = ;

    // Last valid configuration
    // const auto &kd = mat->diffuseColor;
    // return (kd / vec3(M_PI)) * LdotN * lc;
    
    color3 ret = lc * RDM_bsdf(LdotH, NdotH, VdotH, LdotN, VdotN, mat) * LdotN;
    return vec3(min(1.f, ret.x), min(1.f, ret.y), min(1.f, ret.z));
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
    // if (!intersectScene(scene, ray, &intersection))
    // {
    //     return scene->skyColor;
    // }


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

        // bool collision = false;
        // const size_t size = scene->objects.size();
        // for (size_t i = 0; i < size && collision == false; ++i)
        // {
        //     const auto& object = scene->objects[i];
        //     switch (object->geom.type)
        //     {
        //         case PLANE:
        //              collision = (intersectPlane(&shadowRay, &shadowIntersection, object));

        //             break;

        //         case SPHERE:
        //             collision = (intersectSphere(&shadowRay, &shadowIntersection, object));
                
        //             break;

        //         case TRIANGLE:
        //             collision = (intersectTriangle(&shadowRay, &shadowIntersection, object));
                    
        //             break;
            

        //         default : break;
        //     }
        // }

        // if (collision == false)
        // {
        //     color += shade(intersection.normal, -ray->dir, l, light->color, intersection.mat);
        // }
    }


    const auto reflectionDir = reflect(ray->dir, intersection.normal);
    const int add = max(1, 10 - int(intersection.mat->IOR * 10.f));
    
    Ray reflectionRay;
    rayInit(&reflectionRay, intersection.position, reflectionDir, acne_eps, 100000.f, ray->depth + add);

    const color3 cr = trace_ray(scene, &reflectionRay, tree);
    

    const float F = min(1.f, RDM_Fresnel(dot(reflectionRay.dir, intersection.normal), 1.f, intersection.mat->IOR));
    return color + F * cr * intersection.mat->specularColor;
}



color3 superslamping(Scene* scene, KdTree* tree, const vec3& ray_delta_x, const vec3& ray_delta_y, float x, 
                    float y, const vec3& dx, const vec3& dy, int amount) {
    
    color3 result(0);
    for (int xx = 0; xx < amount; ++xx)
    {
        for (int yy = 0; yy < amount; ++yy)
        {
            vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
                (x + (static_cast<float>(xx) / amount)) * dx + (y + (static_cast<float>(yy) / amount)) * dy;

            Ray rx;
            rayInit(&rx, scene->cam.position, normalize(ray_dir));

            result += trace_ray(scene, &rx, tree);
        }
    }

    return result / float(amount * amount);

}

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

float gaussianCoeff(size_t x, size_t y, int matrixSize)
{
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


void computeHeightMap(Scene* scene, KdTree* tree, const vec3& ray_delta_x, 
        const vec3& ray_delta_y,const vec3& dx, const vec3& dy) 
{

    const size_t width  = SceneParameters::imageWidth;
    const size_t height = SceneParameters::imageHeight;
    heightMap.reserve(width * height);

    #pragma omp parallel for
    for (size_t y = 0; y < height; ++y)
    {
        for (size_t x = 0; x < width; ++x)
        {
            vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
                float(x) * dx + float(y) * dy;

            Ray rx;
            rayInit(&rx, scene->cam.position, normalize(ray_dir));
            Intersection intersection;

            // if (intersectScene(scene, &rx, &intersection))
            // {
            //     const float distance = rx.tmax;
            //     heightMap[y * width + x] = min(8.f, abs(SceneParameters::focalDistance - distance));
            // }
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
                for (size_t xx = xxmin; xx < xxmax; ++xx)
                {
                    if (heightMap[yy * img->width + xx] > 0.f)
                    {
                        const float gaussian_coeff = gaussianCoeff((xx - xxmin), (yy - yymin), matrixSize);

                        accumulator += img->data[yy * img->width + xx] * gaussian_coeff;
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

    for (size_t y = 0; y < img->height; ++y)
    {
        for (size_t x = 0; x < img->width; ++x)
        {
            *getPixelPtr(img, x, y) = blurredImage[y * img->width + x];
        }
    }
}



void renderImage(Image *img, Scene *scene)
{
    //! This function is already operational, you might modify it for antialiasing
    //! and kdtree initializaion
    float aspect = 1.f / scene->cam.aspect;

    KdTree *tree = initKdTree(scene);

    //! \todo initialize KdTree

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

#ifdef HEIGHTMAP

    computeHeightMap(scene, tree, ray_delta_x, ray_delta_y, dx, dy);

    const size_t size = img->width * img->height;
    for (size_t i = 0; i < size; ++i)
    {
        color3 *ptr = getPixelPtr(img, i % img->width, i / img->width);
        *ptr = vec3(((heightMap[i] + 8.f) / 2.f) / 8.f);
    }
    return ;
#endif

#ifdef MOTION_BLUR
    std::vector<color3> buff(img->height * img->width, color3(0));
    setCamera(scene, scene->cam.position + vec3(-0.4f, 0, 0), vec3(0,0.3,0), vec3(0, 1, 0), 60, 800.f / 600.f);
    for (size_t a = 0; a < 32; ++a)
    {
        printf("Motion blur, pass: %d\n", a);
#endif 
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

                // Previous code replaced by superslamping
                // vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
                //                              float(i) * dx + float(j) * dy;

                // Ray rx;
                // rayInit(&rx, scene->cam.position, normalize(ray_dir));
                // *ptr = trace_ray(scene, &rx, tree);

                *ptr = superslamping(scene, tree, ray_delta_x, ray_delta_y, i, j, dx, dy, SceneParameters::antiAliasing);
            }

            
#ifdef MOTION_BLUR
        }
        for (int b = 0; b < buff.size(); ++b)
        {
            buff[b] += img->data[b];
        }
        setCamera(scene, scene->cam.position + vec3(0.025f, 0, 0), vec3(0,0.3,0), vec3(0, 1, 0), 60, 800.f / 600.f);
#endif
    }

#ifdef MOTION_BLUR
    for (int b = 0; b < buff.size(); ++b)
    {
        img->data[b] = buff[b] / 32.f;
    }
#endif


    if (SceneParameters::focalDistance > 0)
    {
        printf("Computing depth of field..\n");
        gaussianblur(img, scene, tree, ray_delta_x, ray_delta_y, dx, dy);
    }
}
