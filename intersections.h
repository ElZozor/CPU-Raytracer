#ifndef INTERSECTIONS_H
#define INTERSECTIONS_H

#include "scene_types.h"
#include "scene.h"
#include "ray.h"
#include "raytracer.h"
#include <tuple>


inline color3 getColorFromUV(const Material* mat, const vec2& uv) {
    assert(mat->image != nullptr);
    const auto& image = mat->image;
    const auto& height = mat->height;
    const auto& width  = mat->width;

    const auto uvmod = clamp(uv, 0.0f, 1.0f);
    const auto y = (int((height-1)* uvmod.y+0.5f)) % height;
    const auto x = (int((width-1) * uvmod.x+0.5f)) % width;

    return image[y * width + x];
}


inline bool intersectTriangle(Ray *ray, Intersection *intersection, Object *obj)
{
    const float EPSILON = 1e-7f;
    const auto &triangle = obj->geom.triangle;
    const auto [v0, v1, v2] = triangle.v;

    const vec3 edge0 = v1 - v0;
    const vec3 edge2 = v2 - v0;
    const vec3 h = cross(ray->dir, edge2);

    const float a = dot(edge0, h);

    if (a > -EPSILON && a < EPSILON)
    {
        return false;
    }

    const float f = 1.0f / a;
    const vec3 s = ray->orig - v0;
    const float u = f * dot(s, h);
    if (u < 0.0 || u > 1.0)
    {
        return false;
    }

    const vec3 q = cross(s, edge0);
    const float v = f * dot(ray->dir, q);
    if (v < 0.0 || u + v > 1.0)
    {
        return false;
    }

    float t = f * dot(edge2, q);

    if (t < ray->tmin || t > ray->tmax)
    {
        return false;
    }

    const vec3 N = normalize(cross(edge0, edge2));

    intersection->position = ray->orig + ray->dir * t;
    intersection->normal = N;
    intersection->mat = &obj->mat;
    ray->tmax = t;

    const auto& [width, height] = std::make_tuple(obj->mat.width, obj->mat.height);
    intersection->textured = (width * height != 0);


    if (intersection->textured)
    {
        const auto& [t0, t1, t2] = triangle.vt;
        const float w = 1.f - u - v;
        const auto uv = (t0*w + t1*u + t2*v)/(u + v + w);

        intersection->vt = uv;
    }


    return true;
}








inline bool intersectPlane(Ray *ray, Intersection *intersection, Object *obj)
{
    const auto &plane = obj->geom.plane;
    const auto &d = ray->dir;
    const auto &O = ray->orig;
    const auto &D = plane.dist;
    const auto &n = plane.normal;

    if (dot(d, n) == 0)
    {
        return false;
    }

    const float t = -((dot(O, n) + D) / dot(d, n));

    if (t < ray->tmin || t > ray->tmax)
    {
        return false;
    }

    intersection->mat = &obj->mat;
    intersection->normal = plane.normal;
    intersection->position = (ray->dir * t) + ray->orig;

    ray->tmax = t;

    return true;
}



inline bool intersectSphere(Ray *ray, Intersection *intersection, Object *obj)
{
    const auto &sphere = obj->geom.sphere;
    const auto &R = sphere.radius;
    const auto &d = ray->dir;
    const auto OC = ray->orig - sphere.center;

    const auto b = 2.f * dot(d, OC);
    const auto c = (dot(OC, OC) - R * R);

    const auto delta = (b * b) - 4.f * c;
    if (delta < 0)
    {
        return false;
    }

    float t = -1;
    if (delta == 0)
    {
        t = (-b) / 2.f;
    }
    else
    {
        float t1 = ((-b) - sqrtf(delta)) / 2.f;
        float t2 = ((-b) + sqrtf(delta)) / 2.f;

        if (t1 > t2)
        {
            const float t3 = t1;
            t1 = t2;
            t2 = t3;
        }

        if (t1 > ray->tmin && t1 < ray->tmax)
        {
            t = t1;
        }
        else
        {
            t = t2;
        }
    }

    if (t < ray->tmin || t > ray->tmax)
    {
        return false;
    }

    intersection->mat = &obj->mat;
    intersection->position = (ray->dir * t) + ray->orig;
    intersection->normal = normalize(intersection->position - sphere.center);

    ray->tmax = t;

    return true;
}



inline std::pair<bool, float> intersectCircle(Ray *ray, vec3 center, float radius)
{
    const auto &R = radius;
    const auto &d = ray->dir;
    const auto OC = ray->orig - center;

    const auto b = 2.f * dot(d, OC);
    const auto c = (dot(OC, OC) - R * R);

    const auto delta = (b * b) - 4.f * c;
    if (delta < 0)
    {
        return std::make_pair(false, 0.f);
    }

    float t = -1;
    if (delta == 0)
    {
        t = (-b) / 2.f;
    }
    else
    {
        float t1 = ((-b) - sqrtf(delta)) / 2.f;
        float t2 = ((-b) + sqrtf(delta)) / 2.f;

        if (t1 > t2)
        {
            const float t3 = t1;
            t1 = t2;
            t2 = t3;
        }

        if (t1 > ray->tmin && t1 < ray->tmax)
        {
            t = t1;
        }
        else
        {
            t = t2;
        }
    }

    point3 P = ray->orig + ray->dir * t;
    if (t < ray->tmin || t > ray->tmax || P.y != center.y)
    {
        return std::make_pair(false, 0.f);
    }

    return std::make_pair(true, t);
}




inline bool intersectCone(Ray* ray, Intersection* intersection, Object* object)
{
    const auto& cone = object->geom.cone;
    const auto& pos = ray->orig;
    const auto& dir = ray->dir;
    const auto& center = cone.C;
    const auto& height = cone.V.y;
    const auto& radius = height * std::tan(cone.teta);

    const float A = pos.x - center.x;
    const float B = pos.z - center.z;
    const float D = height - pos.y + center.y;

    const float tan = (radius / height) * (radius / height);

    const float a = (dir.x * dir.x) + (dir.z * dir.z) - (tan*(dir.y * dir.y));
    const float b = (2*A*dir.x) + (2*B*dir.z) + (2*tan*D*dir.y);
    const float c = (A*A) + (B*B) - (tan*(D*D));

    const float delta = b*b - 4*(a*c);
	  if(fabs(delta) < 0.001 || delta < 0.f)
    {
        return false;
    }

    float t1 = (-b - sqrt(delta))/(2*a);
    float t2 = (-b + sqrt(delta))/(2*a);
    float t;

    if (t1 < t2)
    {
        t = t1;
        if (t < ray->tmin || t > ray->tmax)
        {
            return false;
        }
    }
    else
    {
        t = t2;
    }

    const auto [result, value] = intersectCircle(ray, center + cone.V, radius);
    if (result)
    {
        if (value < t && value > ray->tmin && value < ray->tmax)
        {
            t = value;
        }
    }

    if (t > ray->tmax || t < ray->tmin)
    {
        return false;
    }

    const auto P = ray->orig + (ray->dir * t);
    if (P.y > center.y + height || P.y < center.y - height)
    {
        return false;
    } 


    ray->tmax = t;
    intersection->position = P;
    intersection->mat = &object->mat;
    intersection->normal = normalize(intersection->position - center);


     return t;
}


// inline bool intersectCone(Ray* ray, Intersection* intersection, Object* object)
// {
//     const auto& cone = object->geom.cone;
//     const auto& [O, D] = std::make_pair(ray->orig, ray->dir);
//     const auto& C = cone.C, V = cone.V;
//     const auto& teta = cone.teta;

//     const auto DdotV = dot(D,V);
//     const auto cosT = std::cos(teta);
//     const auto costTSquared = cosT * cosT;
//     const auto COdotV = dot(C*O, V);

//     const float a = (DdotV*DdotV) - costTSquared;
//     const float b = 2.f*(DdotV * COdotV - dot(D, C*O*costTSquared));
//     const float c = (COdotV * COdotV) - dot((C*O), (C*O*costTSquared));

//     const float delta = (b*b) - 4.f*a*c;

//     if (delta < 0.f)
//     {
//         return false;
//     }

//     float t;
//     if (delta == 0.f)
//     {
//         t = (-b) / (2.f * a);
//     }
//     else
//     {
//         float t1 = ((-b) -sqrtf(delta)) / (2.f * a);
//         float t2 = ((-b) +sqrtf(delta)) / (2.f * a);

//         if (t1 > t2)
//         {
//             std::swap(t1, t2);
//         }

//         if (t1 > ray->tmin && t1 < ray->tmax)
//         {
//             t = t1;
//         }
//         else
//         {
//             t = t2;
//         }
//     }

//     const vec3 P = ray->orig + (ray->dir * t);
//     if (dot((P-C),V) != cosT)
//     {
//         const float y = dot((P-C),V) - cosT;
//         if (y < 1)
//         printf("%f\n", y);
//         return false;
//     }

//     if (t < ray->tmin || t > ray->tmax)
//     {
//         return false;
//     }
//     printf("okok\n");


//     intersection->mat = &object->mat;
//     intersection->position = ray->orig + (ray->dir * t);

//     ray->tmax = t;

//     return true;
// // }





inline bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection)
{
    bool intersect = false;

    for (const auto &object : scene->objects)
    {
        bool res = false;
        switch (object->geom.type)
        {
        case SPHERE:
            res = intersectSphere(ray, intersection, object);
            break;

        case PLANE:
            res = intersectPlane(ray, intersection, object);
            break;

        case TRIANGLE:
            res = intersectTriangle(ray, intersection, object);
            break;

        case CONE :
            res = intersectCone(ray, intersection, object);
            break;

        default:
            break;
        }

        intersect |= res;
    }

    return intersect;
}


#endif
