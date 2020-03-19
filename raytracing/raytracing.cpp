#include "raytracing.h"
#include "../picture/picture.h"

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <future>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

//TODO : description
/// reflect
/// \param I
/// \param N
/// \return
Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N * 2.f * (I * N);
}

namespace BackGround
{
    int width;
    int height;
    unsigned char *pixmap;
    std::vector<Vec3f> envmap;
};


//TODO : description
/// refract
/// \param I
/// \param N
/// \param eta_t
/// \param eta_i
/// \return
Vec3f refract(const Vec3f &I, const Vec3f &N, const float eta_t, const float eta_i = 1.f) { // Snell's law
    float cosi = -std::max(-1.f, std::min(1.f, I * N));
    if (cosi < 0)
        return refract(I, -N, eta_i, eta_t); // if the ray comes from the inside the object, swap the air and the media
    float eta = eta_i / eta_t;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? Vec3f(1, 0, 0) : I * eta + N * (eta * cosi -
                                                   sqrtf(k)); // k<0 = total reflection, no ray to refract. I refract it anyways, this has no physical meaning
}

bool scene_intersect(const raytracing::entities::Ray &ray,
                     const std::vector<std::unique_ptr<const raytracing::entities::Figure>> &figures,
                     Vec3f &hit, Vec3f &N, raytracing::entities::Material &material) {
    float spheres_dist = std::numeric_limits<float>::max();
    float triangles_dist = std::numeric_limits<float>::max();
    float cubes_dist = std::numeric_limits<float>::max();

    for (const auto &p : figures) {
        float dist_i;
        if (p->ray_intersect(ray, dist_i)) {
            auto &figure_dist = p->NeededDist(spheres_dist, triangles_dist, cubes_dist);
            if (dist_i < figure_dist) {
                figure_dist = dist_i;
                p->SetNeededNormHitMaterial(ray, dist_i, N, hit, material);
            }
        }
    }

    float checkerboard_dist = std::numeric_limits<float>::max();
    if (fabs(ray.dir.y) > 1e-3) {
        float d = -(ray.orig.y + 4) / ray.dir.y; // the checkerboard plane has equation y = -4
        Vec3f pt = ray.orig + ray.dir * d;
        switch (picture::Picture::scene_id) {
            case 1: {
                if (d > 0 && d < spheres_dist && d < triangles_dist) {
                    checkerboard_dist = d;
                    hit = pt;
                    N = Vec3f(0, 1, 0);
                    material.diffuse_color =
                            (int(.5 * hit.x + 1000) + int(.5 * hit.z)) & 1 ? Vec3f(.10, .10, .10) : Vec3f(.3, .2, .1);
                }
            }
            case 2: {
                if (d > 0 && fabs(pt.x) < 10 && pt.z < -10 && pt.z > -30 && d < spheres_dist && d < triangles_dist) {
                    checkerboard_dist = d;
                    hit = pt;
                    N = Vec3f(0, 1, 0);
                    material.diffuse_color =
                            (int(.5 * hit.x + 1000) + int(.5 * hit.z)) & 1 ? Vec3f(.10, .10, .10) : Vec3f(.3, .2, .1);

                }
            }
        }
    }
    return std::min({spheres_dist, checkerboard_dist, cubes_dist, triangles_dist}) < 1000;
}


namespace raytracing {

namespace entities {

Vec3f casting_ray::cast_ray(const Ray &ray,
                            const std::vector<std::unique_ptr<const entities::Figure>> &figures,
                            const std::vector<Light> &lights, size_t depth) {
    Vec3f point, N;
    Material material;

    if (depth > 4 || !scene_intersect(ray, figures, point, N, material)) {
        if (picture::Picture::scene_id==1)
            return Vec3f(0, float(127.0 / 255), float(255.0 / 255));
        else
        {
            int a = std::max(0, std::min(BackGround::width -1, static_cast<int>((atan2(ray.dir.z, ray.dir.x)/(2*M_PI) + .5)*BackGround::width)));
            int b = std::max(0, std::min(BackGround::height-1, static_cast<int>(acos(ray.dir.y)/M_PI*BackGround::height)));
            return BackGround::envmap[a+b*BackGround::width];
        }
    } // background color

    ///TODO make just 2 variables
    Vec3f reflect_dir = reflect(ray.dir, N).normalize();
    Vec3f refract_dir = refract(ray.dir, N, material.refractive_index).normalize();
    Ray reflect_ray(reflect_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3,
                    reflect_dir, &figures, &lights,
                    depth + 1);// offset the original point to avoid occlusion by the object itself (in first param)
    Ray refract_ray(refract_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3,
                    refract_dir, &figures, &lights, depth + 1);

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i = 0; i < lights.size(); ++i) {
        Vec3f light_dir = (lights[i].position - point).normalize();
        float light_distance = (lights[i].position - point).norm();

        Vec3f shadow_orig = light_dir * N < 0 ? point - N * 1e-3 : point + N *
                                                                           1e-3; // checking if the point lies in the shadow of the lights[i]
        Vec3f shadow_pt, shadow_N;
        Material tmp_material;
        if (scene_intersect(Ray(shadow_orig, light_dir), figures,
                            shadow_pt, shadow_N, tmp_material) &&
            (shadow_pt - shadow_orig).norm() < light_distance)
            continue;

        diffuse_light_intensity += lights[i].intensity * std::max(0.f, light_dir * N);
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N) * ray.dir),
                                         material.specular_exponent) * lights[i].intensity;
    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] +
           Vec3f(1., 1., 1.) * specular_light_intensity * material.albedo[1] +
           reflect_ray.color.value() * material.albedo[2] +
           refract_ray.color.value() * material.albedo[3];
}

bool Sphere::ray_intersect(const Ray &ray, float &t0) const {
    Vec3f L = center - ray.orig;
    float tca = L * ray.dir;
    float d2 = L * L - tca * tca;
    if (d2 > radius * radius) return false;
    float thc = sqrtf(radius * radius - d2);
    t0 = tca - thc;
    float t1 = tca + thc;
    if (t0 < 0) t0 = t1;
    if (t0 < 0) return false;
    return true;
}

//TODO needs testing
bool Cube::ray_intersect(const Ray &ray, float &t0) const {
    float t_min, t_max, t_y_min, t_y_max, t_z_min, t_z_max;

    t_min = (bounds[ray.sign[0]].x - ray.orig.x) * ray.invdir.x;
    t_max = (bounds[1 - ray.sign[0]].x - ray.orig.x) * ray.invdir.x;
    t_y_min = (bounds[ray.sign[1]].y - ray.orig.y) * ray.invdir.y;
    t_y_max = (bounds[1 - ray.sign[1]].y - ray.orig.y) * ray.invdir.y;

    if ((t_min > t_y_max) || (t_y_min > t_max))
        return false;
    if (t_y_min > t_min)
        t_min = t_y_min;
    if (t_y_max < t_max)
        t_max = t_y_max;

    t_z_min = (bounds[ray.sign[2]].z - ray.orig.z) * ray.invdir.z;
    t_z_max = (bounds[1 - ray.sign[2]].z - ray.orig.z) * ray.invdir.z;

    if ((t_min > t_z_max) || (t_z_min > t_max))
        return false;
    if (t_z_min > t_min)
        t_min = t_z_min;
    if (t_z_max < t_max)
        t_max = t_z_max;

    t0 = t_min;

    if (t0 < 0) {
        t0 = t_max;
        if (t0 < 0) return false;
    }

    return true;
}

bool Triangle::ray_intersect(const Ray &ray, float &t0) const {
    const float EPSILON = 0.0000001;
    Vec3f vertex0 = p0;
    Vec3f vertex1 = p1;
    Vec3f vertex2 = p2;
    Vec3f edge1, edge2, h, s, q;
    float a, f, u, v;
    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;
    h = cross(ray.dir, edge2);
    a = edge1 * h;
    if (a > -EPSILON && a < EPSILON)
        return false; // This ray is parallel to this triangle.
    f = 1.0 / a;
    s = ray.orig - vertex0;
    u = f * (s * h);
    if (u < 0.0 || u > 1.0)
        return false;
    q = cross(s, edge1);
    v = f * (ray.dir * q);
    if (v < 0.0 || u + v > 1.0)
        return false;
    // At this stage we can compute t to find out where the intersection point is on the line.
    float t = f * (edge2 * q);
    if (t > EPSILON) // ray intersection
    {
        t0 = t;
        return true;
    } else { // This means that there is a line intersection but not a ray intersection.
        return false;
    }
}
}// namespace entities




void render(const char *out_file_path, const std::vector<std::unique_ptr<const entities::Figure>> &figures,
            const std::vector<entities::Light> &lights) {
    const int width = 1920;
    const int height = 1080;
    int env_width,env_height,n=-1;
    const float fov = M_PI / 3.0; ///that's a viewing angle = pi/3
    std::vector<Vec3f> framebuffer(width * height);
    unsigned char *pixmap = stbi_load("../3d.jpg", &env_width, &env_height, &n, STBI_rgb);
    BackGround::width=env_width;
    BackGround::height=env_height;
    BackGround::envmap = std::vector<Vec3f>(BackGround::width*BackGround::height);
    if (!pixmap || 3!=n) {
        std::cerr << "Error: can not load the environment map" << std::endl;
        exit(-1);
    }
    for (int j = BackGround::height-1; j>=0 ; j--) {
        for (int i = 0; i < BackGround::width; i++) {
            BackGround::envmap[i + j * BackGround::width] =
                    Vec3f(pixmap[(i + j * BackGround::width) * 3 + 0], pixmap[(i + j * BackGround::width) * 3 + 1],
                          pixmap[(i + j * BackGround::width) * 3 + 2]) * (1 / 255.);
        }
    }
    stbi_image_free(pixmap);
    const auto amount_of_threads = std::thread::hardware_concurrency(); //because of asynchronius tasks we can make it a bit bigger
//    const auto amount_of_threads = 1; //because of asynchronius tasks we can make it a bit bigger
    std::vector<std::future<void>> tasks(amount_of_threads);
    size_t portion = height / amount_of_threads;
    for (size_t start = 0, finish = portion, index = 0; index < amount_of_threads; ++index, start = finish, finish = index == amount_of_threads - 1 ? height : finish + portion) {
        tasks[index] = std::async(std::launch::async, [&, start, finish, width, height, fov](){
           for (size_t j = start; j < finish; ++j){
               for (size_t i = 0; i < width; ++i){
                   auto dir_x = (i + 0.5) - width / 2.;
                   auto dir_y = -(j + 0.5) + height / 2.;    // this flips the image at the same time
                   auto dir_z = -height / (2. * tan(fov / 2.));
                   framebuffer[i + j * width] = entities::casting_ray::cast_ray(entities::Ray(
                            Vec3f(0, 0, 0), Vec3f(dir_x, dir_y, dir_z).normalize()), figures, lights);
               }
           }
        });
    }
    for (auto &&task : tasks) {
        task.get();
    }

//    std::ofstream ofs; // save the framebuffer to file
//    ofs.open(out_file_path, std::ios::binary);
//    ofs << "P6\n" << width << " " << height << "\n255\n";
    std::vector<unsigned char> Npixmap(width*height*3);
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c*(1./max);
        for (size_t j = 0; j<3; j++) {
            Npixmap[i*3+j] = (unsigned char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    stbi_write_jpg("FINAL.jpg", width, height, 3, Npixmap.data(), 100);
}

}// namespace raytracing