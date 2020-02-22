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
    const auto& t0 = triangle.t0;
    const auto& t1 = triangle.t1;
    const auto& t2 = triangle.t2;
    const float minx = min(t0.x, min(t1.x, t2.x));
    const float maxx = max(t0.x, max(t1.x, t2.x));
    const float miny = min(t0.y, min(t1.y, t2.y));
    const float maxy = max(t0.y, max(t1.y, t2.y));
    const float minz = min(t0.z, min(t1.z, t2.z));
    const float maxz = max(t0.z, max(t1.z, t2.z));

    return {
        vec3(minx, miny, minz),
        vec3(maxx, maxy, maxz)
    };
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
    bool intersectsRay(Ray* ray)
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
    bool intersectsTriangle(const Object* object)
    {
        const AABB tb = triangleBounds(object);
        return !((bounds.min.x > tb.max.x)      // trop à droite
                || (bounds.max.x < tb.min.x)    // trop à gauche
                || (bounds.min.y > tb.max.y)    // trop en bas
                || (bounds.max.y < tb.min.y)    // trop en haut
                || (bounds.min.z > tb.max.z)    // derrière
                || (bounds.max.z < tb.min.z));  // devant
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
                if (intersectsTriangle(object))
                {
                    result = true;
                }
                else
                {
                    result = false;
                }
                break;

            case SPHERE:
            {
                const auto& sphere = object->geom.sphere;
                if (intersectsSphere(sphere.center, sphere.radius))
                {
                    result = true;
                }
                else
                {
                    result = false;
                }
                break;
            }

            default :
                printf("Unhandled object insertion \n");
                exit(1);
                break; 
        }

        return result;
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

                        default: break;
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
                        case SPHERE:
                            result = intersectSphere(ray, intersection, object) || result;
                            break;

                        case PLANE :
                            result = intersectPlane(ray, intersection, object) || result;
                            break;

                        case TRIANGLE :
                            result = intersectTriangle(ray, intersection, object) || result;
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