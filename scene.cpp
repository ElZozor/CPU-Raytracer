#include "scene.h"
#include "scene_types.h"
#include <string.h>
#include <algorithm>

Object *initSphere(point3 center, float radius, Material mat) {
    Object *ret;
    ret = (Object *)malloc(sizeof(Object));
    ret->geom.type = SPHERE;
    ret->geom.sphere.center = center;
    ret->geom.sphere.radius = radius;
    memcpy(&(ret->mat), &mat, sizeof(Material));
    return ret;
}

Object *initPlane(vec3 normal, float d, Material mat) {
    Object *ret;
    ret = (Object *)malloc(sizeof(Object));
    ret->geom.type = PLANE;
    ret->geom.plane.normal = normalize(normal);
    ret->geom.plane.dist = d;
    memcpy(&(ret->mat), &mat, sizeof(Material));
    return ret;
}
#include "stdio.h"
Object *initTriangle(point3 v[3], vec3 vt[3], vec3 vn[3], Material mat) {
    Object *ret;
    
    ret = (Object *)malloc(sizeof(Object));
    ret->geom.type = TRIANGLE;

    memcpy(ret->geom.triangle.v,  v,  sizeof(point3) * 3);
    memcpy(ret->geom.triangle.vt, vt, sizeof(vec3)   * 3);
    memcpy(ret->geom.triangle.vn, vn, sizeof(vec3)   * 3);
    memcpy(&(ret->mat), &mat, sizeof(Material));

    assert(v[0].x == ret->geom.triangle.v[0].x);
    printf("%f %f %f \n",v[0].x, v[1].x, v[2].x);
    // printf("%f %f %f \n",ret->geom.triangle.v[0], ret->geom.triangle.v[1], ret->geom.triangle.v[2]);

    return ret;
}


Object *initCone(point3 C, vec3 V, float teta, Material mat) {
    Object *ret;
    ret = (Object *)malloc(sizeof(Object));
    ret->geom.type = CONE;
    ret->geom.cone.C = C;
    ret->geom.cone.V = V;
    ret->geom.cone.teta = teta;

    memcpy(&(ret->mat), &mat, sizeof(Material));

    return ret;
}


void freeObject(Object *obj) {
    free(obj);
}

Light *initLight(point3 position, color3 color) {
    Light *light = (Light*)malloc(sizeof(Light));
    light->position = position;
    light->color = color;
    return light;
}

void freeLight(Light *light) {
    free(light);
}

Scene * initScene() {
    return new Scene;
}

void freeScene(Scene *scene) {
    std::for_each(scene->objects.begin(), scene->objects.end(), freeObject);
    std::for_each(scene->lights.begin(), scene->lights.end(), freeLight);
    delete scene;
}

void setCamera(Scene *scene, point3 position, point3 at, vec3 up, float fov, float aspect) {
    scene->cam.fov = fov;
    scene->cam.aspect = aspect;
    scene->cam.position = position;
    scene->cam.zdir = normalize(at-position);
    scene->cam.xdir = normalize(cross(up, scene->cam.zdir));
    scene->cam.ydir = normalize(cross(scene->cam.zdir, scene->cam.xdir));
    scene->cam.center = 1.f / tanf ((scene->cam.fov * glm::pi<float>() / 180.f) * 0.5f) * scene->cam.zdir;
}

void addObject(Scene *scene, Object *obj) {
    scene->objects.push_back(obj);
}

void addLight(Scene *scene, Light *light) {
    scene->lights.push_back(light);
}

void setSkyColor(Scene *scene, color3 c) {
    scene->skyColor = c;
}
