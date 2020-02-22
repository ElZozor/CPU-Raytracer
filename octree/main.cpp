#include "octree.hpp"

int main(void)
{
    Octree octree(true);

    Material mat;
    Object* sphere = initSphere(point3(0,0,0), 10.f, mat);

    return 0;
}