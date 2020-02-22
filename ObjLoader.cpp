#include "ObjLoader.hpp"
#include <fstream>
#include <iostream>
#include <sstream>

ObjLoader::ObjLoader()
{

}

void ObjLoader::load(const std::string& filename)
{
    Material mat = {1.1022, 0.0579, {0.286, 0.235, 0.128}, {1.0, 0.766, 0.762}};
    std::vector<point3> points;

    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Unable to open file : " << filename << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(file, line))
    {
        if (line.size() == 0 || line.find(" ") == std::string::npos)
        {
            continue;
        }


        std::istringstream iss(line);
        std::string id;
        iss >> id;
        if (id == "#")
        {
            continue;
        }
        else if (id == "mtllib")
        {
            std::string filename;
            iss >> filename;
            loadMaterialFile(filename);
        }
        else if (id == "v")
        {
            float x,y,z,w = 1.f;
            iss >> x >> y >> z >> w;

            points.push_back(point3(x,y,z));
        }
        else if (id == "usemtl")
        {
            std::string name;
            iss >> name;
            if (mMaterials.find(name) != mMaterials.end())
            {
                mat = mMaterials[name];
            }
        }
        else if (id == "f")
        {
            int garbage;
            int vs[4];
            for (size_t i = 0; i < 3; ++i)
            {
                iss >> vs[i]; iss.ignore(2); iss >> garbage;
                if (vs[i] < 0)
                {
                    vs[i] = points.size() + vs[i];
                }
                
                --vs[i];
            }

            mObjects.push_back(
                initTriangle(
                    points[vs[0]],
                    points[vs[1]],
                    points[vs[2]],
                    mat
                )
            );
        }
    }



    file.close();
}


void ObjLoader::loadMaterialFile(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cout << "Unable to find file: " << filename << std::endl;
        return;
    }

    Material mat;
    bool constructingMaterial = false;
    std::string materialName;

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string id;
        iss >> id;


        if (id == "newmtl")
        {
            if (constructingMaterial)
            {
                mMaterials[materialName] = mat;
            }
            
            iss >> materialName;
            constructingMaterial = true;
        } 
        else if (id == "Kd")
        {
            float a,b,c;
            iss >> a >> b >>c;
            mat.diffuseColor = color3(a,b,c);
        }
        else if (id == "Ks")
        {
            float a,b,c;
            iss >> a >> b >> c;
            mat.specularColor = color3(a,b,c);
        }
        else if (id == "Ni")
        {
            float a;
            iss >> a;
            mat.roughness = a;
        }
        else if (id == "Ka")
        {
            float r;
            iss >> r;
            mat.IOR = r;
        }
    }

    if (constructingMaterial)
    {
        mMaterials[materialName] = mat;
    }


    file.close();
}


std::vector<Object*>& ObjLoader::getObjects()
{
    return mObjects;
}