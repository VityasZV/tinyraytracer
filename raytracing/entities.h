//
// Created by Виктор Задябин on 09.03.2020.
//

#ifndef TINYRAYTRACER_ENTITIES_H
#define TINYRAYTRACER_ENTITIES_H

#include "../geometry.h"
#include <optional>
#include <array>
#include <memory>
namespace raytracing::entities {
enum class FigureName{
    Sphere,
    Triangle,
    Cube
};

class Figure; class Sphere; struct Ray; struct Light; class Cube; class Triangle;
namespace casting_ray {

    //TODO : description
    /// cast_ray - casts new rays for reflect, refract.
    /// \param ray
    /// \param spheres
    /// \param lights
    /// \param depth
    /// \return
    Vec3f cast_ray(const Ray &ray, const std::vector<std::unique_ptr<const Figure>> &figures,
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
    void SetColor(const std::vector<std::unique_ptr<const Figure>> &figures, const std::vector<Light> &lights, size_t depth){
        color = casting_ray::cast_ray(*this, figures, lights, depth);
    }

public:
    Vec3f orig, dir;       // ray orig and dir
    Vec3f invdir;          // inverse direction
    std::array<int, 3> sign;
    std::optional<Vec3f> color = std::nullopt;


    /// Ray
    /// \param orig - source of ray
    /// \param dir - direction of ray
    Ray(const Vec3f &orig, const Vec3f &dir, const std::vector<std::unique_ptr<const Figure>> * const figures_ptr = nullptr,
            const std::vector<Light>* const lights_ptr = nullptr,
            const std::optional<size_t> depth = std::nullopt) : orig(orig), dir(dir) {
        invdir = 1 / dir;
        sign[0] = (invdir.x < 0);
        sign[1] = (invdir.y < 0);
        sign[2] = (invdir.z < 0);
        if (depth.has_value()) {
            SetColor(*figures_ptr, *lights_ptr, depth.value());
        }
    }
    Ray(const Ray& r) : orig(r.orig), dir(r.dir), invdir(r.invdir), sign(r.sign), color(r.color){}
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

class Figure {
private:
    Material material;
    FigureName name;
public:
    Figure(const Material &m, const FigureName &n): material(m), name(n){}
    Figure(const Figure &f): material(f.material), name(f.name) {}
    virtual ~Figure() = default;
    virtual bool ray_intersect(const Ray& ray, float &t0) const = 0;
    virtual float& NeededDist(const float &spheres_dist, const float &triangles_dist, const float &cubes_dist) const = 0;
    virtual void SetNeededNormHitMaterial(const Ray& ray, const float& dist_i, Vec3f &n, Vec3f &h, Material &m) const = 0;
    Material GetMaterial() const {return material;}
    FigureName GetFigureName() const {return name;}
};

class Triangle : public Figure {
private:
    Vec3f p0, p1, p2;
public:
    Triangle(Vec3f p0, Vec3f p1, Vec3f p2, const Material &m) : Figure(m, FigureName::Triangle), p0(p0), p1(p1), p2(p2) {}

    Triangle(const std::vector<Vec3f> &vec, const Material &m) : Figure(m, FigureName::Triangle), p0(vec[0]), p1(vec[1]), p2(vec[2]) {}


    bool ray_intersect(const Ray &ray, float &t0) const override;

    float& NeededDist(const float &spheres_dist, const float &triangles_dist, const float &cubes_dist) const override {
        return const_cast<float &>(triangles_dist);
    }

    //нормаль
    Vec3f get_N() const{
        return (cross(p1 - p0, p2 - p0)).normalize();
    }
    void SetNeededNormHitMaterial(const Ray& ray, const float& dist_i, Vec3f &n, Vec3f &h, Material &m) const override {
        n = this->get_N();
        h = ray.orig + ray.dir * dist_i;
        m = GetMaterial();
    }
};


class Sphere : public Figure {
private:
    Vec3f center;
    float radius;
public:
    /// Sphere - object that constructed by
    /// \param c - its center
    /// \param r - its radius
    /// \param m - its material
    Sphere(const Vec3f &c, const float r, const Material &m) : Figure(m, FigureName::Sphere), center(c), radius(r) {}
    ~Sphere() = default;
    Sphere(const Sphere &s) : Figure(s.GetMaterial(), s.GetFigureName()), center(s.center), radius(s.radius){}
    /// ray_intersect - checks if ray intersects sphere
    /// \param r - Ray
    /// \param t0 - coordinate of possible intersection
    /// \return
    bool ray_intersect(const Ray& r, float &t0) const override;
    float& NeededDist(const float &spheres_dist, const float &triangles_dist, const float &cubes_dist) const override {
        return const_cast<float &>(spheres_dist);
    }
    void SetNeededNormHitMaterial(const Ray& ray, const float& dist_i, Vec3f &n, Vec3f &h, Material &m) const override {
        h = ray.orig + ray.dir * dist_i;
        n = (h - center).normalize();
        m = GetMaterial();
    }
};

class Cube : public Figure {
private:
    ///The bounds of the volume define a set of lines parallel to each axis
    ///of the coordinate system which we can also expressed using the line equation.
    std::array<Vec3f, 2> bounds;
public:
    ///To represent an axis-aligned bounding volume, all we need are two points
    /// representing the minimum and maximum extent of the box (called bounds in the code).
    Cube(const Vec3f &vmin, const Vec3f &vmax, const Material &m) : Figure(m, FigureName::Cube){
        bounds[0] = vmin;//левая нижняя ближняя?
        bounds[1] = vmax;//правая верхняя дальная?
    }
    ~Cube() = default;
    Cube(const Cube &c) : Figure(c.GetMaterial(), c.GetFigureName()), bounds(c.bounds){}
    bool ray_intersect(const Ray &ray, float &t0) const override;
    float& NeededDist(const float &spheres_dist, const float &triangles_dist, const float &cubes_dist) const override {
        return const_cast<float &>(cubes_dist);
    }
    void SetNeededNormHitMaterial(const Ray& ray, const float& dist_i, Vec3f &n, Vec3f &h, Material &m) const override {
        h = ray.orig + ray.dir * dist_i;
        //TODO - неверно, надо разобраться
        n = (h - bounds[0]).normalize();
        m = GetMaterial();
    }
};

}
#endif //TINYRAYTRACER_ENTITIES_H
