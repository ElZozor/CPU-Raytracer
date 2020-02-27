#include "kdtree.h"
#include "scene_types.h"
#include "intersections.h"


#include <stdio.h>
#include <vector>
#include <algorithm>
#include <functional>


#define MAX_OBJECT 32
#define MAX_RECURSION 256


/**
 * @brief Compute the minimum values for each axis of the vector
 * 
 * @param a         The first vector
 * @param b         The second vector
 * @return vec3     The minimum for each fields of them
 */
inline vec3 vec3min(const vec3& a, const vec3& b)
{
    return vec3(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z));
}

/**
 * @brief Compute the maximum values for each axis of the vector
 * 
 * @param a         The first vector
 * @param b         The second vector
 * @return vec3     The maximum for each fields of them
 */
inline vec3 vec3max(const vec3&a, const vec3& b)
{
    return vec3(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z));
}


/**
 * @brief Compute the sphere AABB
 * 
 * @param object    The sphere
 * @return AABB     The AABB coordinates
 */
inline static AABB sphereBounds(const Object* object)
{
    const auto& sphere = object->geom.sphere;

    return {
        sphere.center - vec3(sphere.radius),
        sphere.center + vec3(sphere.radius)
    };
}


/**
 * @brief Compute the triangle AABB bounds
 * 
 * @param object    The triangle
 * @return AABB     The AABB coordinates
 */
inline static AABB triangleBounds(const Object* object)
{
    const auto& triangle = object->geom.triangle;
    const auto& [v0, v1, v2] = triangle.v;
    const float minx = min(v0.x, min(v1.x, v2.x));
    const float maxx = max(v0.x, max(v1.x, v2.x));
    const float miny = min(v0.y, min(v1.y, v2.y));
    const float maxy = max(v0.y, max(v1.y, v2.y));
    const float minz = min(v0.z, min(v1.z, v2.z));
    const float maxz = max(v0.z, max(v1.z, v2.z));

    return {
        vec3(minx, miny, minz),
        vec3(maxx, maxy, maxz)
    };
}


inline static AABB coneBounds(const Object* object)
{
    const auto& cone = object->geom.cone;
    const auto& height = cone.V.y;
    const auto& radius = height * std::tan(cone.teta);

    printf("height: %f\n", radius);

    const float miny = cone.C.y - height;
    const float maxy = cone.C.y + height;
    const float minx = cone.C.x - 2.f*radius;
    const float maxx = cone.C.x + 2.f*radius;
    const float minz = cone.C.z - 2.f*radius;
    const float maxz = cone.C.z + 2.f*radius;



    return {vec3(minx, miny, minz), vec3(maxx, maxy, maxz)};
}


/**
 * @brief Return the object bounds 
 * 
 * @param object 
 * @return AABB 
 */
static AABB objectBounds(const Object* object)
{
    assert(object->geom.type != PLANE);

    switch (object->geom.type)
    {
        case TRIANGLE :
            return triangleBounds(object);

        case SPHERE:
            return sphereBounds(object);

        case CONE:
            return coneBounds(object);

        default:
            printf("Object geometry not handled !\n");
            exit(1);
    }
}


struct KdTreeNode
{
    KdTreeNode* left;
    KdTreeNode* right;

    AABB bounds;
    std::vector<Object*> objects;

    std::function<std::pair<float, float>(const Object*, int)> composantGetter = {
        [](const Object* a, int dim) { 
            const auto bounds = objectBounds(a);
            return std::pair<float, float>(bounds.min[dim], bounds.max[dim]);
        }
    };

    bool isLeaf;
    float splitPoint;

    static int nodenumber;
    int number;
    int dim;
    

    KdTreeNode(AABB cbounds, int cdim) :
        left(nullptr),
        right(nullptr),
        bounds(cbounds),
        isLeaf(true),
        number(nodenumber++),
        dim(cdim)
    {
        objects.reserve(MAX_OBJECT);
    }


    /**
     * @brief Compute the collision between current node and AABB
     * 
     * @param aabb
     * @return true 
     * @return false 
     */
    inline bool intersectsAABB(const AABB& aabb)
    {
        return !((bounds.min.x > aabb.max.x)      // trop à droite
                || (bounds.max.x < aabb.min.x)    // trop à gauche
                || (bounds.min.y > aabb.max.y)    // trop en bas
                || (bounds.max.y < aabb.min.y)    // trop en haut
                || (bounds.min.z > aabb.max.z)    // derrière
                || (bounds.max.z < aabb.min.z));  // devant
    }


    /**
     * @brief Return if the current node intersects a sphere or not
     * 
     * @param sphereCenter  The sphere center
     * @param sphereRadius  The sphere radius
     * @return true 
     * @return false 
     */
    bool intersectsSphere(const vec3& sphereCenter, float sphereRadius)
    {
        const auto& [mini, maxi] = std::make_tuple(bounds.min, bounds.max);
        vec3 closestPointInAabb = min(max(sphereCenter, mini), maxi);
        vec3 seg = closestPointInAabb - sphereCenter;
        float distanceSquared = dot(seg, seg);
        // The AABB and the sphere overlap if the closest point within the rectangle is
        // within the sphere's radius
        return distanceSquared < (sphereRadius * sphereRadius);
    }

    /**
     * @brief Return whether or not the current node intersects this ray
     * 
     * @param ray       The ray to test
     * @return true 
     * @return false 
     */
    inline bool intersectsRay(Ray* ray)
    {
        const auto& [mini, maxi, orig, dir] = std::make_tuple(bounds.min, bounds.max, ray->orig, ray->dir);
        float t[8];
        t[0] = (mini.x - orig.x) / dir.x;
        t[1] = (maxi.x - orig.x) / dir.x;
        t[2] = (mini.y - orig.y) / dir.y;
        t[3] = (maxi.y - orig.y) / dir.y;
        t[4] = (mini.z - orig.z) / dir.z;
        t[5] = (maxi.z - orig.z) / dir.z;
        t[6] = fmax(fmax(fmin(t[0], t[1]), fmin(t[2], t[3])), fmin(t[4], t[5]));
        t[7] = fmin(fmin(fmax(t[0], t[1]), fmax(t[2], t[3])), fmax(t[4], t[5]));

        return (t[7] >= 0) && (t[6] <= t[7]);
    }


    /**
     * @brief Return whether or not the current node intersects this triangle
     * 
     * @param object    The triangle
     * @return true 
     * @return false 
     */
    inline bool intersectsTriangle(const Object* object)
    {
        return intersectsAABB(triangleBounds(object));
    }


    /**
     * @brief Return whether or not the current node intersects this cone
     * 
     * @param object    The cone
     * @return true 
     * @return false 
     */
    inline bool intersectsCone(const Object* object)
    {
        return intersectsAABB(coneBounds(object));
    }



    /**
     * @brief Return whether or not the current node intersects the given object
     * 
     * @param object    The object to test
     * @return true 
     * @return false 
     */
    bool intersectsObject(const Object* object)
    {
        bool result = false;
      
        switch (object->geom.type)
        {
            case TRIANGLE : 
                result = intersectsTriangle(object);
                break;

            case SPHERE:
            {
                const auto& sphere = object->geom.sphere;
                result = intersectsSphere(sphere.center, sphere.radius);
                break;
            }

            case CONE:
                result = intersectsCone(object);
                break;

            default :
                printf("Unhandled object insertion \n");
                exit(1);
                break; 
        }

        return result;
    }



    inline std::tuple<float, float, float> computeRayTs(Ray* ray) 
    {
        const auto& [dmin, dmax, dO, dD] = std::make_tuple(
            bounds.min[dim], bounds.max[dim], ray->orig[dim], ray->dir[dim]
        );
        
        return std::make_tuple(
            (dmin - dO) / dD,
            (splitPoint - dO) / dD,
            (dmax - dO) / dD
        );
    }


    /**
     * @brief Compute the nearest ray intersection
     * 
     * @param ray           The ray
     * @param intersection  The intersection which will be updated
     * @return true 
     * @return false 
     */
    bool computeRayIntersections(Ray *ray, Intersection *intersection)
    {
        if (intersectsRay(ray))
        {
            if (isLeaf)
            {
                bool result = false;
                for (const auto& object : objects)
                {
                    switch (object->geom.type)
                    {
                        case SPHERE:
                            result = intersectSphere(ray, intersection, object) || result;
                            break;

                        case PLANE :
                            result = intersectPlane(ray, intersection, object) || result;
                            break;

                        case TRIANGLE :
                            result = intersectTriangle(ray, intersection, object) || result;
                            break;
                        
                        case CONE:
                            result = intersectCone(ray, intersection, object) || result;
                            break;

                        default: 
                            printf("Unhandled ray/object intersection !\n");
                            exit(3);
                    }
                }

                return result;
            }
            else
            {
                // assert(left->intersectsRay(ray) || right->intersectsRay(ray));
                

                bool result = left->computeRayIntersections(ray, intersection);
                result = right->computeRayIntersections(ray, intersection) || result;
                
                

                return result;
            }
        }

        return false;
    }



    /**
     * @brief Return whether or not the ray intersects an object \
     * Stop on the first collision 
     * 
     * @param ray               The ray
     * @param intersection      The intersection which will be updated
     * @return true 
     * @return false 
     */
    bool doesRayIntersectsAnObject(Ray* ray, Intersection* intersection)
    {
        if (intersectsRay(ray))
        {
            if (isLeaf)
            {
                bool result = false;
                for (const auto& object : objects)
                {
                    switch (object->geom.type)
                    {
                        case SPHERE :
                            result = intersectSphere(ray, intersection, object) || result;
                            break;

                        case PLANE :
                            result = intersectPlane(ray, intersection, object) || result;
                            break;

                        case TRIANGLE :
                            result = intersectTriangle(ray, intersection, object) || result;
                            break;
                        
                        case CONE :
                            result = intersectCone(ray, intersection, object) || result;
                            break;

                        default: break;
                    }

                    if (result)
                    {
                        return true;
                    }
                }
            }
            else
            {
                // assert(left->intersectsRay(ray) || right->intersectsRay(ray));
                return left->computeRayIntersections(ray, intersection)
                    || right->computeRayIntersections(ray, intersection);
            }
        }

        return false;
    }



    /**
     * @brief Intsert an object if it fit on the current node
     * 
     * @param object    The object to insert
     * @return true 
     * @return false 
     */
    bool insert(Object* object)
    {
        bool result;

        if (isLeaf)
        {
            if (intersectsObject(object))
            {
                objects.push_back(object);
                
                if (objects.size() == MAX_OBJECT && number < MAX_RECURSION)
                {
                    split();
                }
            }



            return true;
        }
        else
        {
            result = left ->insert(object);
            result = right->insert(object) || result;
        }

        return result;
    }


    /**
     * @brief Split the current node on it's current axis
     * 
     */
    void split()
    {
        if (dim == 0)
        {
            splitXAxis();
        } 
        else if (dim == 1)
        {
            splitYAxis();
        }
        else
        {
            splitZAxis();
        }

        isLeaf = false;
        for (const auto& object: objects)
        {
            const bool insertLeft  = left->insert(object);
            const bool insertRight = right->insert(object);
            assert(insertLeft || insertRight);
        }

        objects = std::vector<Object*>(); //Faster than vector.clear()
    }


    /**
     * @brief Sort the current objects given the current node dimension (X,Y,Z)
     * 
     */
    void sortObjects()
    {
        std::sort(objects.begin(), objects.end(), 
            [&composantGetter = composantGetter, &dim = dim](const Object* a, const Object* b) {
                return composantGetter(a, dim).first < composantGetter(b, dim).first;
            }
        );
    }


    /**
     * @brief Split the node on the X axis
     * 
     */
    void splitXAxis() 
    {
        sortObjects();
        splitPoint = composantGetter(objects[MAX_OBJECT / 2], dim).first;
        // splitPoint = (bounds.max.x + bounds.min.x)/2.f;

        left  = new KdTreeNode({bounds.min, vec3(splitPoint, bounds.max.y, bounds.max.z)}, (dim + 1) % 3);
        right = new KdTreeNode({vec3(splitPoint, bounds.min.y, bounds.min.z), bounds.max}, (dim + 1) % 3);
    }

    /**
     * @brief Split the node on the Y axis
     * 
     */
    void splitYAxis() 
    {
        sortObjects();
        splitPoint = composantGetter(objects[MAX_OBJECT / 2], dim).first;
        // splitPoint = (bounds.max.y + bounds.min.y)/2.f;

        left  = new KdTreeNode({bounds.min, vec3(bounds.max.x, splitPoint, bounds.max.z)}, (dim + 1) % 3);
        right = new KdTreeNode({vec3(bounds.min.x, splitPoint, bounds.min.z), bounds.max}, (dim + 1) % 3);
    }

    /**
     * @brief Split the node on the Z axis
     * 
     */
    void splitZAxis() 
    {
        sortObjects();
        splitPoint = composantGetter(objects[MAX_OBJECT / 2], dim).first;
        // splitPoint = (bounds.max.z + bounds.min.z)/2.f;

        left  = new KdTreeNode({bounds.min, vec3(bounds.max.x, bounds.max.y, splitPoint)}, (dim + 1) % 3);
        right = new KdTreeNode({vec3(bounds.min.x, bounds.min.y, splitPoint), bounds.max}, (dim + 1) % 3);
    }

    
};

int KdTreeNode::nodenumber = 0;




struct s_kdtree
{
    KdTreeNode *root;
    std::vector<Object*> outOfTree;

    /**
     * @brief Construct a new s_kdtree object
     * 
     * @param scene     The scene
     */
    s_kdtree(const Scene* scene)
    {
        root = new KdTreeNode(computeSceneBounds(scene), 0);
        insertObjects(scene);
    }


    /**
     * @brief Compute the bounds of the current scene
     * 
     * @param scene     The scene
     * @return AABB     The bounds of the scene
     */
    AABB computeSceneBounds(const Scene* scene) const
    {
        AABB treeBounds;
        bool first = true;
        for (const auto& object : scene->objects)
        {
            AABB bounds;

            switch (object->geom.type)
            {
                case PLANE:
                    continue;
                    break;

                case SPHERE:
                    bounds = sphereBounds(object);
                    break;
                
                case TRIANGLE: 
                    bounds = triangleBounds(object);
                    break;
                
                case CONE :
                    bounds = coneBounds(object);
                    break;

                default:
                    printf("KD-TREE: unknow object type !\n %d", object->geom.type);
                    exit(1);
                    break;
            }

            if (first)
            {
                treeBounds = bounds;
                first = false;
            }
            else
            {
                treeBounds.min = vec3min(treeBounds.min, bounds.min);
                treeBounds.max = vec3max(treeBounds.max, bounds.max);
            } 
        }

        return treeBounds;
    }



    /**
     * @brief Insert all object into the scene
     * 
     * @param scene     The scene
     */
    void insertObjects(const Scene* scene)
    {
        for (const auto& object: scene->objects)
        {
            switch (object->geom.type)
            {
                case PLANE :
                    outOfTree.push_back(object);
                    break;
                
                case SPHERE:
                    root->insert(object);
                    break;
                
                case TRIANGLE:
                    root->insert(object);
                    break;
                
                case CONE:
                    root->insert(object);
                    break;
                
                default:
                    printf("Unhandled object insertion\n");
                    exit(2);
            }
        }
    }


    /**
     * @brief Compute the nearest object/ray intersection
     * 
     * @param ray           The ray
     * @param intersection  The output intersection
     * @return true 
     * @return false 
     */
    bool computeIntersection(Ray *ray, Intersection *intersection)
    {
        bool intersect = false;
        for (const auto& object : outOfTree)
        {
            intersect = intersectPlane(ray, intersection, object) || intersect;
        }

        return root->computeRayIntersections(ray, intersection) || intersect;
    }


    /**
     * @brief Return whether or not the ray intersects an object
     * 
     * @param ray               The ray
     * @param intersection      The output intersection
     * @return true 
     * @return false 
     */
    bool intersectAnObject(Ray* ray, Intersection* intersection)
    {
        for (const auto& object : outOfTree)
        {
            if (intersectPlane(ray, intersection, object))
            {
                return true;
            }
        }

        return root->doesRayIntersectsAnObject(ray, intersection);
    }
};



/**
 * @brief Simply call the kdtree constructor
 * Insert all the objects into it and return the resulting tree
 * 
 * @param scene     The scene
 * @return KdTree*  The resulting tree
 */
KdTree* initKdTree(Scene *scene)
{
    return new KdTree(scene);
}



/**
 * @brief Return the nearest object intersection
 * 
 * @param scene             The current scene
 * @param tree              The concerned kdtree
 * @param ray               The ray
 * @param intersection      The output intersection
 * @return true
 * @return false 
 */
bool intersectKdTree(Scene *scene, KdTree *tree, Ray *ray, Intersection *intersection)
{
    return tree->computeIntersection(ray, intersection);
}



/**
 * @brief Return whether or not the ray intersects an object
 * 
 * @param scene             The current scene
 * @param tree              The concerned kdtree
 * @param ray               The ray
 * @param intersection      The output intersection
 * @return true 
 * @return false 
 */
bool intersectObjectKdTree(Scene *scene, KdTree *tree, Ray *ray, Intersection *intersection)
{
    return tree->intersectAnObject(ray, intersection);
}