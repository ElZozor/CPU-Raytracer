#ifndef OCTREE_HPP
#define OCTREE_HPP

#define CAPACITY 8

#include <set>
#include <array>
#include "../scene_types.h"

using namespace glm;

struct ObjectComparator {
    bool operator() (const Object* a, const Object* b) const {
        return a < b;
    }
};

class Octree
{
private : //Variables
    bool root;
    bool leaf;

    vec3 position;
    vec3 size;

    std::set<const Object*, ObjectComparator> objects;
    Octree *childs;
public:
    Octree(bool root);
    ~Octree();

    bool insert(const Object* object);

private:
    Octree();
    void split();
    void setBounds(std::pair<vec3, vec3> bounds);
    // void updateBounds(const Object* object);
    bool intersects(const Object* object) const;
    bool intersectsSphere(const vec3 center, const float radius) const;
    std::array<std::pair<vec3, vec3>, 8> computeChildsBounds();

};

#endif // OCTREE_HPP
