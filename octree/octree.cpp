#include "octree.hpp"

Octree::Octree(bool croot) :
    root(croot),
    leaf(true)
{

}

Octree::Octree() :
    root(false),
    leaf(true)
{

}

Octree::~Octree()
{
    if (!leaf)
    {
        delete [] childs;
    }
}

bool Octree::insert(const Object *object)
{
    if (!intersects(object))
    {
        return false;

        // updateBounds(object);
    }

    bool inserted(false);
    if (leaf)
    {
        if (objects.insert(object).second)
        {
            if (objects.size() == CAPACITY)
            {
                split();
            }

            inserted = true;
        }
        else
        {
            inserted = false;
        }
    }
    else
    {
        for (size_t i = 0; i < 8U; ++i)
        {
            inserted = childs[i].insert(object) || inserted;
        }
    }

    return inserted;
}


void Octree::split()
{
    leaf = false;
    childs = new Octree[8]();

    const auto bounds = computeChildsBounds();
    for (size_t i = 0; i < 8U; ++i)
    {
        childs[i].setBounds(bounds[i]);
        
        for (const auto object: objects)
        {
            childs[i].insert(object);
        }
    }
    
    objects.clear();
}

void Octree::setBounds(const std::pair<vec3, vec3> bounds)
{
    position = bounds.first;
    size     = bounds.second;
}

bool Octree::intersects(const Object *object) const
{
    switch(object->geom.type)
    {
        case SPHERE :
        {
            return intersectsSphere(object->geom.sphere.center, object->geom.sphere.radius);
        }
            

        case PLANE :
            return false;

        default :
            return false;
    }
}

bool Octree::intersectsSphere(const vec3 center, const float radius) const
{
    const vec3 maximum = (position + size);
    const float x = max(position.x, min(center.x, maximum.x));
    const float y = max(position.y, min(center.y, maximum.y));
    const float z = max(position.z, min(center.y, maximum.z));

    const float distance = sqrt((x - center.x) * (x - center.x) +
                            (y - center.y) * (y - center.y) +
                            (z - center.z) * (z - center.z));

    return distance < radius;
}



std::array<std::pair<vec3, vec3>, 8> Octree::computeChildsBounds()
{
    const vec3 updatedSize = size / 2.f;

    return std::array<std::pair<vec3, vec3>, 8>({
                std::pair<vec3, vec3>(position + vec3(0,             0,             0),             updatedSize),
                std::pair<vec3, vec3>(position + vec3(updatedSize.x, 0,             0),             updatedSize),
                std::pair<vec3, vec3>(position + vec3(0,             updatedSize.y, 0),             updatedSize),
                std::pair<vec3, vec3>(position + vec3(updatedSize.x, updatedSize.y, 0),             updatedSize),
                std::pair<vec3, vec3>(position + vec3(0,             0,             updatedSize.z), updatedSize),
                std::pair<vec3, vec3>(position + vec3(updatedSize.x, 0,             updatedSize.z), updatedSize),
                std::pair<vec3, vec3>(position + vec3(0,             updatedSize.y, updatedSize.z), updatedSize),
                std::pair<vec3, vec3>(position + vec3(updatedSize.x, updatedSize.y, updatedSize.z), updatedSize)
    });
}





















