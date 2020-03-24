#include "raytracing.h"
#include "entities.h"


#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <future>
#include <memory>

void tree_trace(std::shared_ptr<raytracing::kd_tree::KdTree::Node> tree) {
    if (tree == nullptr) return;
    std::cout << "\nOBHOD\n" << std::endl;
    std::cout << tree->box.Size() << tree->box.GetVMin() << tree->box.GetVMax() << std::endl;
    auto *obj = std::get_if<const std::vector<std::shared_ptr<raytracing::kd_tree::KdTree::RenderWrapper>>>(
            &tree->plane_or_figures);
    if (obj) {
        std::cout << "всего фигур: " << obj->size() << " штук: " << "vmin=" << tree->box.GetVMin() << "vmax= "
                  << tree->box.GetVMax() << std::endl;
        for (auto &p : *obj) {
            p->obj->print();
        }
    } else {
        auto *pl = std::get_if<raytracing::entities::Plane>(&tree->plane_or_figures);
        std::cout << pl->GetPos();
    }
    tree_trace(tree->child.first);
    tree_trace(tree->child.second);
}

std::unordered_map<std::shared_ptr<const raytracing::entities::Figure>, float>
raytracing::Render::bin_search_in_tree(const raytracing::entities::Ray &ray,
                                       std::shared_ptr<raytracing::kd_tree::KdTree::Node> tree) {
    float t_near, t_far;
    if (tree == nullptr) return {}; //if tree is not formed at all
    if (tree->box.Intersect(ray, t_near, t_far)) {
        auto *objects = std::get_if<const std::vector<std::shared_ptr<raytracing::kd_tree::KdTree::RenderWrapper>>>(
                &tree->plane_or_figures);
        if (objects) {
            std::unordered_map<std::shared_ptr<const raytracing::entities::Figure>, float> result;
            for (const auto &p : *objects) {
                float dist_i;
                if (p->obj->ray_intersect(ray, dist_i)) {
                    result.insert({p->obj, dist_i});
                }
            }
            return result;
        }
        auto res1 = tree->child.first ? bin_search_in_tree(ray, tree->child.first)
                                      : std::unordered_map<std::shared_ptr<const raytracing::entities::Figure>, float>();
        auto res2 = tree->child.second ? bin_search_in_tree(ray, tree->child.second)
                                       : std::unordered_map<std::shared_ptr<const raytracing::entities::Figure>, float>();
        auto result = [&res1, &res2]() {
            std::unordered_map<std::shared_ptr<const raytracing::entities::Figure>, float> res;
            for (auto &e : res1){
                res.insert(e);
            }
            for (auto &e : res2){
                res.insert(e);
            }
            return res;
        }();
        return result;

        //returning result
        //FIXME upper code goes throughout all nodes of a tree, checking only those nodes,
        // which have aabb intersection with ray, and if we find figures there->we are checking
        // figures on intersection
        //FIXME now it works a bit incorrectly because a dont know which figure intersected first -> wrong dist of sphere/triangle -> bad shadows and wrong
        // show of figures -- FIXED with vector of nodes that ray intersects. right now i deleted down code -> maybe not finally
    }
    return {};
}

void raytracing::Render::initialize_kd_tree(std::shared_ptr<raytracing::kd_tree::KdTree::Node> tree) {
    raytracing::tree = tree;
}

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
    float checkerboard_dist = std::numeric_limits<float>::max();

    //if we launched tree
    if (raytracing::tree) {
        auto needed_figures_and_dists = raytracing::Render::bin_search_in_tree(ray, raytracing::tree);
        if (needed_figures_and_dists.size() != 0) {
            for (auto &f : needed_figures_and_dists){
                auto &figure_dist = f.first->NeededDist(spheres_dist, triangles_dist, cubes_dist);
                if (f.second < figure_dist) {
                    figure_dist = f.second;
                    f.first->SetNeededNormHitMaterial(ray, f.second, N, hit, material);
                }
            }
        } else {
            //ray doesnt intersect with any of primitives
            //need to find checkerboard
            if (fabs(ray.dir.y) > 1e-3) {
                float d = -(ray.orig.y + 4) / ray.dir.y; // the checkerboard plane has equation y = -4
                Vec3f pt = ray.orig + ray.dir * d;
                if (d > 0 && d < spheres_dist && d < triangles_dist && d < cubes_dist) {
                    checkerboard_dist = d;
                    hit = pt;
                    N = Vec3f(0, 1, 0);
                    material.diffuse_color =
                            (int(.5 * hit.x + 1000) + int(.5 * hit.z)) & 1 ? Vec3f(.10, .10, .10) : Vec3f(.3, .2, .1);
                }
            }
        }
    } else {
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
        if (fabs(ray.dir.y) > 1e-3) {
            float d = -(ray.orig.y + 4) / ray.dir.y; // the checkerboard plane has equation y = -4
            Vec3f pt = ray.orig + ray.dir * d;
            if (d > 0 && d < spheres_dist && d < triangles_dist && d < cubes_dist) {
                checkerboard_dist = d;
                hit = pt;
                N = Vec3f(0, 1, 0);
                material.diffuse_color =
                        (int(.5 * hit.x + 1000) + int(.5 * hit.z)) & 1 ? Vec3f(.10, .10, .10) : Vec3f(.3, .2, .1);
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
        return Vec3f(0, float(127.0 / 255), float(255.0 / 255)); // background color
    }

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
    float t1 = (bounds[0].x - ray.orig.x)*ray.invdir.x;
    float t2 = (bounds[1].x - ray.orig.x)*ray.invdir.x;
    float t3 = (bounds[0].y - ray.orig.y)*ray.invdir.y;
    float t4 = (bounds[1].y - ray.orig.y)*ray.invdir.y;
    float t5 = (bounds[0].z - ray.orig.z)*ray.invdir.z;
    float t6 = (bounds[1].z - ray.orig.z)*ray.invdir.z;
    float t = 0;
    float t_near = std::max(std::max(std::min(t1, t2), std::min(t3, t4)), std::min(t5, t6));
    float t_far = std::min(std::min(std::max(t1, t2), std::max(t3, t4)), std::max(t5, t6));

    // if tmax < 0, ray (line) is intersecting AABB, but the whole AABB is behind us
    if (t_far < 0)
    {
        t = t_far;
        return false;
    }

    // if tmin > tmax, ray doesn't intersect AABB
    if (t_near> t_far)
    {
        t = t_far;
        return false;
    }

    t = t_near;
    //to is nearest point of intersection
    t0 = t_near;
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


Vec3f anti_aliasing(double dir_x, double dir_y, double dir_z,
                    const std::vector<std::unique_ptr<const entities::Figure>> &figures,
                    const std::vector<entities::Light> &lights) {
    Vec3f anti_alias = Vec3f(0, 0, 0);
    for (int k = 0; k < 5; ++k) {
        switch (k % 5) {
            case 0:
                anti_alias = anti_alias + entities::casting_ray::cast_ray(entities::Ray(
                        Vec3f(0, 0, 0), Vec3f(dir_x, dir_y, dir_z).normalize()), figures, lights);
                break;
            case 1:
                anti_alias = anti_alias + entities::casting_ray::cast_ray(entities::Ray(
                        Vec3f(0, 0, 0), Vec3f(dir_x + 0.5, dir_y, dir_z).normalize()), figures, lights);
                break;
            case 2:
                anti_alias = anti_alias + entities::casting_ray::cast_ray(entities::Ray(
                        Vec3f(0, 0, 0), Vec3f(dir_x, dir_y + 0.5, dir_z).normalize()), figures, lights);
                break;
            case 3:
                anti_alias = anti_alias + entities::casting_ray::cast_ray(entities::Ray(
                        Vec3f(0, 0, 0), Vec3f(dir_x, dir_y, dir_z + 0.5).normalize()), figures, lights);
                break;
            default:
                anti_alias = anti_alias + entities::casting_ray::cast_ray(entities::Ray(
                        Vec3f(0, 0, 0), Vec3f(dir_x, dir_y, dir_z - 0.5).normalize()), figures, lights);
                break;
        }
    }
    return (anti_alias / 5);
}


void Render::render(const char *out_file_path, const std::vector<std::unique_ptr<const entities::Figure>> &figures,
                    const std::vector<entities::Light> &lights) {
    //tree_trace(raytracing::tree);
    const int width = 1920;
    const int height = 1080;
    const float fov = M_PI / 3.0; ///that's a viewing angle = pi/3
    std::vector<Vec3f> framebuffer(width * height);
    const auto amount_of_threads = std::thread::hardware_concurrency(); //because of asynchronius tasks we can make it a bit bigger
    std::vector<std::future<void>> tasks(amount_of_threads);
    size_t portion = height / amount_of_threads;
    for (size_t start = 0, finish = portion, index = 0; index < amount_of_threads; ++index, start = finish, finish =
            index == amount_of_threads - 1 ? height : finish + portion) {
        tasks[index] = std::async(std::launch::async, [&, start, finish, width, height, fov]() {
            for (size_t j = start; j < finish; ++j) {
                for (size_t i = 0; i < width; ++i) {
                    auto dir_x = (i + 0.5) - width / 2.;
                    auto dir_y = -(j + 0.5) + height / 2.;    // this flips the image at the same time
                    auto dir_z = -height / (2. * tan(fov / 2.));
                    //FIXME needs check on scene 3, if scene == 3 then without anti_aliasing
                    framebuffer[i + j * width] = entities::casting_ray::cast_ray(entities::Ray(
                            Vec3f(0, 0, 0), Vec3f(dir_x, dir_y, dir_z).normalize()), figures, lights);
                    //framebuffer[i + j * width] = raytracing::anti_aliasing(dir_x, dir_y, dir_z, figures, lights);
                }
            }
        });
    }
    for (auto &&task : tasks) {
        task.get();
    }

    std::ofstream ofs; // save the framebuffer to file
    ofs.open(out_file_path, std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height * width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max > 1) c = c * (1. / max);
        for (size_t j = 0; j < 3; ++j) {
            ofs << (char) (255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}


}// namespace raytracing
