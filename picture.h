//
// Created by Виктор Задябин on 10.03.2020.
//

#ifndef RAYTRACER_PICTURE_H
#define RAYTRACER_PICTURE_H

#include "raytracing/entities.h"

#include <unordered_map> 

namespace picture {

namespace {

enum class MaterialName {
    ivory,
    glass,
    mirror,
    red_rubber
};
struct SphereParams {
    SphereParams(const Vec3f& coordinates, const float &radius, const raytracing::entities::Material &material) :
            coordinates(coordinates), material(material), radius(radius){}
    Vec3f coordinates;
    raytracing::entities::Material material;
    float radius;
};

struct LightParams {
    LightParams(const Vec3f &position, const float &intensity) :
            position(position), intensity(intensity){}
    Vec3f position;
    float intensity;
};

struct CubeParams {
    CubeParams(const Vec3f &vmin, const Vec3f &vmax, const raytracing::entities::Material &m) :
                vmin(vmin), vmax(vmax), material(m){}
    Vec3f vmin;
    Vec3f vmax;
    raytracing::entities::Material material;
};

}// namespace

class Picture {
private:
    void PreparingOutFileAndScene (int argc, const char** argv){
        std::unordered_map<std::string, std::string> cmd_line_params;
        for(int i=0; i<argc; i++)
        {
            std::string key(argv[i]);
            if(key.size() > 0 && key[0]=='-')
            {
                if(i != argc-1) // not last argument
                {
                    cmd_line_params[key] = argv[i+1];
                    i++;
                }
                else
                    cmd_line_params[key] = "";
            }
        }   
        if(cmd_line_params.find("-out") != cmd_line_params.end())
            out_file_path = cmd_line_params["-out"];
    } 

    std::unordered_map<MaterialName, raytracing::entities::Material> Materials {
            {MaterialName::ivory, raytracing::entities::Material(1.0, Vec4f(0.6, 0.3, 0.1, 0.0), Vec3f(0.4, 0.4, 0.3), 50.)},
            {MaterialName::glass, raytracing::entities::Material(1.5, Vec4f(0.0, 0.5, 0.1, 0.8), Vec3f(0.6, 0.7, 0.8), 125.)},
            {MaterialName::red_rubber, raytracing::entities::Material(1.0, Vec4f(0.9, 0.1, 0.0, 0.0), Vec3f(0.3, 0.1, 0.1), 10.)},
            {MaterialName::mirror, raytracing::entities::Material(1.0, Vec4f(0.0, 10.0, 0.8, 0.0), Vec3f(1.0, 1.0, 1.0), 1425.)}
    };
    std::vector<SphereParams> spheres_params{
            {Vec3f(-3, 0, -16), 2, Materials[MaterialName::ivory]},
            {Vec3f(-1.0, -1.5, -12), 2, Materials[MaterialName::glass]},
            {Vec3f(1.5, -0.5, -18), 3, Materials[MaterialName::red_rubber]},
            {Vec3f(7, 5, -18), 4, Materials[MaterialName::mirror]}
    };

    std::vector<LightParams> lights_params{
            {Vec3f(-20, 20, 20), 1.5},
            {Vec3f(30, 50, -25), 1.8},
            {Vec3f(30, 20, 30), 1.7}
    };

    std::vector<CubeParams> cubes_params{
            {Vec3f(-3, 0, -3), Vec3f(1, 4, -7), Materials[MaterialName::red_rubber]}
    };
public:
    std::string out_file_path = "out.ppm";
    Picture(const int argc, const char** argv){
        PreparingOutFileAndScene (argc, argv);
        for (const auto& p : spheres_params){
            spheres.emplace_back(raytracing::entities::Sphere(p.coordinates, p.radius, p.material));
        }
        for (const auto& p : lights_params){
            lights.emplace_back(raytracing::entities::Light(p.position, p.intensity));
        }
//        for (const auto& p : cubes_params){
//            cubes.emplace_back(raytracing::entities::Cube(p.vmin, p.vmax, p.material));
//        }
    }
    ~Picture() = default;
    std::vector<raytracing::entities::Sphere> spheres;
    std::vector<raytracing::entities::Light> lights;
    std::vector<raytracing::entities::Cube> cubes;
};

}// namespace picture
#endif //RAYTRACER_PICTURE_H
