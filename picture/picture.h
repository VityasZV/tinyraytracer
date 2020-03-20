//
// Created by Виктор Задябин on 10.03.2020.
//

#ifndef RAYTRACER_PICTURE_H
#define RAYTRACER_PICTURE_H

#include "../raytracing/entities.h"

#include <unordered_map>
#include <memory>

namespace picture {

namespace {

enum class MaterialName {
    ivory,
    glass,
    mirror,
    red_rubber
};

struct SphereParams {
    SphereParams(const Vec3f &coordinates, const float &radius, const raytracing::entities::Material &material) :
            coordinates(coordinates), material(material), radius(radius) {}

    Vec3f coordinates;
    raytracing::entities::Material material;
    float radius;
};

struct LightParams {
    LightParams(const Vec3f &position, const float &intensity) :
            position(position), intensity(intensity) {}

    Vec3f position;
    float intensity;
};

struct CubeParams {
    CubeParams(const Vec3f &vmin, const Vec3f &vmax, const raytracing::entities::Material &m) :
            vmin(vmin), vmax(vmax), material(m) {}

    Vec3f vmin;
    Vec3f vmax;
    raytracing::entities::Material material;
};

struct TriangleParams {
    TriangleParams(const Vec3f &p0, const Vec3f &p1, const Vec3f &p2, const raytracing::entities::Material &m) :
            material(m), p0(p0), p1(p1), p2(p2) {}

    Vec3f p0, p1, p2;
    raytracing::entities::Material material;
};

}// namespace

class Picture {
private:
    void PreparingOutFileAndScene(int argc, const char **argv);

    std::unordered_map<MaterialName, raytracing::entities::Material> Materials{
            {MaterialName::ivory,      raytracing::entities::Material(1.0, Vec4f(0.6, 0.3, 0.1, 0.0),
                                                                      Vec3f(0.4, 0.4, 0.3), 50.)},
            {MaterialName::glass,      raytracing::entities::Material(1.5, Vec4f(0.0, 0.5, 0.1, 0.8),
                                                                      Vec3f(0.6, 0.7, 0.8), 125.)},
            {MaterialName::red_rubber, raytracing::entities::Material(1.0, Vec4f(0.9, 0.1, 0.0, 0.0),
                                                                      Vec3f(0.3, 0.1, 0.1), 10.)},
            {MaterialName::mirror,     raytracing::entities::Material(1.0, Vec4f(0.0, 10.0, 0.8, 0.0),
                                                                      Vec3f(1.0, 1.0, 1.0), 1425.)}
    };
    std::vector<SphereParams> spheres_params{
            {Vec3f(-3, 0, -16),      2, Materials[MaterialName::ivory]},
            {Vec3f(-1.0, -1.5, -12), 2, Materials[MaterialName::glass]},
            {Vec3f(1.5, -0.5, -18),  3, Materials[MaterialName::red_rubber]},
            {Vec3f(7, 5, -18),       4, Materials[MaterialName::mirror]},
            {Vec3f(-8, 5, -18),      4, Materials[MaterialName::mirror]}

    };

    std::vector<LightParams> lights_params{
            {Vec3f(-20, 20, 20), 1.5},
            {Vec3f(30, 50, -25), 1.8},
            {Vec3f(30, 20, 30),  1.7}
    };

    std::vector<CubeParams> cubes_params{
            {Vec3f(-3, 0, -3), Vec3f(1, 4, -7), Materials[MaterialName::red_rubber]}
    };
    std::vector<TriangleParams> triangle_params{
            //dont forget about right trio while adding params!!!
            //front
            {Vec3f(-6, 0, -6), Vec3f(-5, 0, -6), Vec3f(-5, 2, -6), Materials[MaterialName::red_rubber]},
            {Vec3f(-5, 2, -6), Vec3f(-6, 2, -6), Vec3f(-6, 0, -6), Materials[MaterialName::red_rubber]},
            //down
            {Vec3f(-5, 0, -6), Vec3f(-6, 0, -6), Vec3f(-5, 0, -8), Materials[MaterialName::red_rubber]},
            {Vec3f(-5, 0, -8), Vec3f(-6, 0, -6), Vec3f(-6, 0, -8), Materials[MaterialName::red_rubber]},
            //up
            {Vec3f(-6, 2, -6), Vec3f(-5, 2, -6), Vec3f(-5, 2, -8), Materials[MaterialName::red_rubber]},
            {Vec3f(-5, 2, -8), Vec3f(-6, 2, -8), Vec3f(-6, 2, -6), Materials[MaterialName::red_rubber]},
            //left
            {Vec3f(-6, 2, -8), Vec3f(-6, 2, -6), Vec3f(-6, 0, -6), Materials[MaterialName::red_rubber]},
            {Vec3f(-6, 0, -8), Vec3f(-6, 2, -8), Vec3f(-6, 0, -6), Materials[MaterialName::red_rubber]},
            //right
            {Vec3f(-5, 0, -8), Vec3f(-5, 2, -8), Vec3f(-5, 2, -6), Materials[MaterialName::red_rubber]},
            {Vec3f(-5, 0, -8), Vec3f(-5, 2, -6), Vec3f(-5, 0, -6), Materials[MaterialName::red_rubber]},

    };

    void MakeTriangleMash(const char *file_name);

public:
    std::string out_file_path = "out.ppm";
    int scene_id = 1;

    Picture(const int argc, const char **argv);

    ~Picture() = default;

    std::vector<raytracing::entities::Light> lights;
    std::vector<std::unique_ptr<const raytracing::entities::Figure>> figures;
};

}// namespace picture
#endif //RAYTRACER_PICTURE_H
