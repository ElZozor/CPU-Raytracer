#include "ObjLoader.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include "SplittableString.hpp"
#include "color_msg.h"

ObjLoader::ObjLoader()
{

}

void ObjLoader::load(const std::string& filename)
{
    mFilename = filename;
    
    Material mat = {1.1022, 0.0579, {0.286, 0.235, 0.128}, {1.0, 0.766, 0.762}};
    std::vector<point3> v({{0,0,0}});
    std::vector<vec3> vt({{0,0,0}}), vn({{0,0,0}});

    std::ifstream file(filename);
    if (!file.is_open())
    {
        logError(0, "Unable to open/find the file");
    }

    SplittableString line;
    std::istringstream iss;
    for (size_t line_number = 0; std::getline(file, line); ++line_number)
    {
        if (line.size() == 0 || line.find(" ") == std::string::npos)
        {
            continue;
        }

        const std::vector<SplittableString> values = line.split(" ");
        if (values.size() < 1)
        {
            continue;
        }

        const auto& id = values[0];
        if (id == "#")
        {
            continue;
        }
        else if (id == "mtllib")
        {
            if (values.size() < 2)
            {
                logError(line_number, "'mtllib', missing filename");
            }
            
            loadMaterialFile(values[1]);
        }
        else if (id == "v")
        {
            if (values.size() < 4)
            {
                logError(line_number, "'v', missing operand");
            }
            
            v.push_back(point3(
                std::stof(values[1]),
                std::stof(values[2]),
                std::stof(values[3])
            ));
        }
        else if (id == "vt")
        {
            vec3 vti(0);
            for (size_t i = 1; i < values.size() && i < 4; ++i)
            {
                vti[i-1] = std::stof(values[i]);
            }

            vt.push_back(vti);
        }
        else if (id == "vn")
        {
            vec3 vni(0);
            for (size_t i = 1; i < values.size() && i < 3; ++i)
            {
                vni[i-1] = std::stof(values[i]);
            }

            vn.push_back(vni);
        }
        else if (id == "usemtl")
        {
            if (values.size() < 2)
            {
                logError(line_number, "'usemtl', missing material name");
            }

            std::string name = values[1];
            if (mMaterials.find(name) != mMaterials.end())
            {
                mat = mMaterials[name];
            }
        }
        else if (id == "f")
        {
            point3 vs[3];
            vec3 vts[3], vns[3];
            for (size_t i = 1; i < 4; ++i)
            {
                auto coordinates = values[i].split("/");
                if (coordinates.size() < 3)
                {
                    logError(line_number, "'f' missing operand");
                }
                

                const auto indice = i-1;
                if (!coordinates[0].empty())
                {
                    int ind = std::stoi(coordinates[0]);
                    if (ind < 0) ind = v.size() - ind;
                    vs[indice]  = v[ind];
                }
                if (!coordinates[1].empty())
                {
                    int ind = std::stoi(coordinates[1]);
                    if (ind < 0) ind = vt.size() - ind;
                    vts[indice] = vt[ind];
                }

                if (!coordinates[2].empty())
                {
                    int ind = std::stoi(coordinates[2]);
                    if (ind < 0) ind = vn.size() - ind;
                    vns[indice] = vn[ind];
                }
            }
            mObjects.push_back(
                initTriangle(
                    vs, vts, vns, mat
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


void ObjLoader::logError(size_t line, const std::string& message)
{
    logc(LIGHT_RED, stderr, "Error in file %s at line %ld: %s\n", mFilename.c_str(), line, message.c_str());
    logc(LIGHT_RED, stderr, "Aborting..\n");
    exit(1);
}