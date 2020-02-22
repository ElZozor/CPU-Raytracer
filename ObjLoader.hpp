#ifndef OBJ_LOADER_HPP
#define OBJ_LOADER_HPP

#include <string>
#include <map>
#include "scene_types.h"

class ObjLoader
{
private:
    std::vector<Object*> mObjects;
    std::map<std::string, Material> mMaterials;

public:
    explicit ObjLoader();
    ~ObjLoader() = default;

    void load(const std::string& filename);
    std::vector<Object*>& getObjects();

private:
    void loadMaterialFile(const std::string& filename);
};


#endif //OBJ_LOADER_HPP