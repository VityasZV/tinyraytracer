//
// Created by Виктор Задябин on 09.03.2020.
//

#ifndef TINYRAYTRACER_ENTITIES_H
#define TINYRAYTRACER_ENTITIES_H

#include "optional"
#include "../geometry.h"


namespace raytracing::entities {

struct Sphere; struct Ray; struct Light; struct Cube;
namespace casting_ray {

    //TODO : description
    /// cast_ray - casts new rays for reflect, refract.
    /// \param ray
    /// \param spheres
    /// \param lights
    /// \param depth
    /// \return
    Vec3f cast_ray(const Ray &ray, const std::vector<Sphere> &spheres,
                   const std::vector<Light> &lights, size_t depth = 0);
}

struct Light {
    Vec3f position;
    float intensity;
    /// Light is a source of light
    /// \param p - position of light source
    /// \param i - intensity of light source
    Light(const Vec3f &p, const float i) : position(p), intensity(i) {}
};

struct Ray {
private:
    /// sets color for ray, calling cast_ray function
    /// \param spheres
    /// \param lights
    /// \param depth
    void SetColor(const std::vector<Sphere> &spheres, const std::vector<Light> &lights, size_t depth){
        color = casting_ray::cast_ray(*this, spheres, lights, depth);
    }

public:
    Vec3f orig, dir;       // ray orig and dir
    Vec3f invdir;          // inverse direction
    int sign[3];
    std::optional<Vec3f> color = std::nullopt;


    /// Ray
    /// \param orig - source of ray
    /// \param dir - direction of ray
    Ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> * const spheres_ptr = nullptr,
            const std::vector<Light>* const lights_ptr = nullptr, const std::optional<size_t> depth = std::nullopt) : orig(orig), dir(dir) {
        invdir = 1 / dir;
        sign[0] = (invdir.x < 0);
        sign[1] = (invdir.y < 0);
        sign[2] = (invdir.z < 0);
        if (depth.has_value()) {
            SetColor(*spheres_ptr, *lights_ptr, depth.value());
        }
    }
};

struct Material {
    float refractive_index;
    Vec4f albedo;
    Vec3f diffuse_color;
    float specular_exponent;

    /// Material is a
    /// \param r
    /// \param a
    /// \param color
    /// \param spec
    Material(const float r, const Vec4f &a, const Vec3f &color, const float spec) : refractive_index(r),
                                                                                    albedo(a),
                                                                                    diffuse_color(color),
                                                                                    specular_exponent(spec) {}

    Material() : refractive_index(1), albedo(1, 0, 0, 0), diffuse_color(), specular_exponent() {}
};

namespace {
    struct Figure {
        Material material;
        Figure(const Material& m): material(m){}
        virtual ~Figure() = default;
        virtual bool ray_intersect(const Ray& ray, float &t0) const = 0;
    };

}// namespace



struct Sphere : public Figure {
    Vec3f center;
    float radius;
    /// Sphere - object that constructed by
    /// \param c - its center
    /// \param r - its radius
    /// \param m - its material
    Sphere(const Vec3f &c, const float r, const Material &m) : Figure(m), center(c), radius(r) {}
    ~Sphere() = default;

    /// ray_intersect - checks if ray intersects sphere
    /// \param r - Ray
    /// \param t0 - coordinate of possible intersection
    /// \return
    bool ray_intersect(const Ray& r, float &t0) const override;
};

struct Cube : public Figure {
    Cube(const Vec3f &vmin, const Vec3f &vmax, const Material &m) : Figure(m){
        bounds[0] = vmin;
        bounds[1] = vmax;
    }
    ~Cube() = default;

    bool ray_intersect(const Ray &ray, float &t0) const override;
    Vec3f bounds[2];
};

}
#endif //TINYRAYTRACER_ENTITIES_H
