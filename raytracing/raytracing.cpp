#include "raytracing.h"

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <future>

//TODO : description
/// reflect
/// \param I
/// \param N
/// \return
Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N * 2.f * (I * N);
}

//TODO : description
/// refract
/// \param I
/// \param N
/// \param eta_t
/// \param eta_i
/// \return
Vec3f refract(const Vec3f &I, const Vec3f &N, const float eta_t, const float eta_i=1.f) { // Snell's law
    float cosi = - std::max(-1.f, std::min(1.f, I*N));
    if (cosi<0) return refract(I, -N, eta_i, eta_t); // if the ray comes from the inside the object, swap the air and the media
    float eta = eta_i / eta_t;
    float k = 1 - eta*eta*(1 - cosi*cosi);
    return k<0 ? Vec3f(1,0,0) : I*eta + N*(eta*cosi - sqrtf(k)); // k<0 = total reflection, no ray to refract. I refract it anyways, this has no physical meaning
}

bool scene_intersect(const raytracing::entities::Ray& ray, const std::vector<raytracing::entities::Sphere> &spheres, Vec3f &hit, Vec3f &N, raytracing::entities::Material &material) {
    float spheres_dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(ray, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = ray.orig + ray.dir*dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }

    float checkerboard_dist = std::numeric_limits<float>::max();
    if (fabs(ray.dir.y)>1e-3)  {
        float d = -(ray.orig.y+4)/ray.dir.y; // the checkerboard plane has equation y = -4
        Vec3f pt = ray.orig + ray.dir*d;
        if (d>0 && d<spheres_dist) {
            checkerboard_dist = d;
            hit = pt;
            N = Vec3f(0,1,0);
            material.diffuse_color = (int(.5*hit.x+1000) + int(.5*hit.z)) & 1 ? Vec3f(.10, .10, .10) : Vec3f(.3, .2, .1);
        }
    }
    return std::min(spheres_dist, checkerboard_dist)<1000;
}


namespace raytracing {

namespace entities{

    Vec3f casting_ray::cast_ray(const Ray &ray,
                                const std::vector<Sphere> &spheres,
                                const std::vector<Light> &lights, size_t depth) {
        Vec3f point, N;
        Material material;

        if (depth>4 || !scene_intersect(ray, spheres, point, N, material)) {
            return Vec3f(0, float(127.0/255), float(255.0/255)); // background color
        }

        ///TODO make just 2 variables
        Vec3f reflect_dir = reflect(ray.dir, N).normalize();
        Vec3f refract_dir = refract(ray.dir, N, material.refractive_index).normalize();
        Ray reflect_ray(reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3,
                                              reflect_dir, &spheres, &lights, depth+1);// offset the original point to avoid occlusion by the object itself (in first param)
        Ray refract_ray(refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3,
                                              refract_dir, &spheres, &lights, depth + 1);

        float diffuse_light_intensity = 0, specular_light_intensity = 0;
        for (size_t i=0; i<lights.size(); i++) {
            Vec3f light_dir      = (lights[i].position - point).normalize();
            float light_distance = (lights[i].position - point).norm();

            Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // checking if the point lies in the shadow of the lights[i]
            Vec3f shadow_pt, shadow_N;
            Material tmp_material;
            if (scene_intersect(Ray(shadow_orig, light_dir), spheres,
                                shadow_pt, shadow_N, tmp_material) &&
                (shadow_pt - shadow_orig).norm() < light_distance)
                continue;

            diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir*N);
            specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*ray.dir),
                                             material.specular_exponent) * lights[i].intensity;
        }
        return material.diffuse_color * diffuse_light_intensity * material.albedo[0] +
               Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1] +
               reflect_ray.color.value()*material.albedo[2] +
               refract_ray.color.value()*material.albedo[3];
    }

    bool Sphere::ray_intersect(const Ray& ray, float &t0) const {
        Vec3f L = center - ray.orig;
        float tca = L*ray.dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }

}// namespace entities




void render(const std::vector<entities::Sphere> &spheres, const std::vector<entities::Light> &lights) {
    const int   width    = 1920;
    const int   height   = 1080;
    const float fov      = M_PI/3.; ///that's a viewing angle = pi/3
    std::vector<Vec3f> framebuffer(width*height);
    const auto amount_of_threads = std::thread::hardware_concurrency();
    std::vector<std::future<void>> tasks(amount_of_threads);

    size_t portion = width / amount_of_threads;
    for (size_t j = 0; j<height; j++) { // actual rendering loop
        //trying to parallel compute pieces of "line" of a picture
        for (size_t start = 0, finish = portion, index = 0;
        index < amount_of_threads; ++index, start=finish,
        finish = index == amount_of_threads - 1 ? width : finish + portion){
            tasks[index] = std::async(std::launch::async, [=, &spheres, &lights, &framebuffer](){
                for (size_t i = start; i < finish; i++) {
                    auto dir_x = (i + 0.5) - width / 2.;
                    auto dir_y = -(j + 0.5) + height / 2.;    // this flips the image at the same time
                    auto dir_z = -height / (2. * tan(fov / 2.));
                    framebuffer[i + j * width] = entities::casting_ray::cast_ray(entities::Ray(
                            Vec3f(0, 0, 0), Vec3f(dir_x, dir_y, dir_z).normalize()),
                                                          spheres, lights);
                }
            });
        }
        /// waiting our parallel tasks
        for (auto&& task : tasks){
            task.get();
        }
    }

    std::ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm",std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c*(1./max);
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

}// namespace raytracing
