#ifndef INTERSECTIONS_H
#define INTERSECTIONS_H

#include "scene_types.h"
#include "scene.h"
#include "ray.h"
#include "raytracer.h"



inline bool intersectTriangle(Ray *ray, Intersection *intersection, Object *obj)
{

    const float EPSILON = 0.0000001f;
    const auto &triangle = obj->geom.triangle;
    const vec3 v0 = triangle.t0;
    const vec3 v1 = triangle.t1;
    const vec3 v2 = triangle.t2;

    const vec3 edge1 = v1 - v0;
    const vec3 edge2 = v2 - v0;
    const vec3 h = cross(ray->dir, edge2);

    const float a = dot(edge1, h);

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

    const vec3 q = cross(s, edge1);
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

    const vec3 N = normalize(cross(edge1, edge2));

    intersection->position = ray->orig + ray->dir * t;
    intersection->normal = N;
    intersection->mat = &obj->mat;

    ray->tmax = t;
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

        default:
            break;
        }

        intersect |= res;
    }

    return intersect;
}


#endif