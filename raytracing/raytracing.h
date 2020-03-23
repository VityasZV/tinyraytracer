#ifndef TINYRAYTRACER_RAYTRACING_H
#define TINYRAYTRACER_RAYTRACING_H

#include "../geometry/geometry.h"
#include "entities.h"
#include <memory>
#include "../kd-tree/kdtree.h"

namespace raytracing {
static std::shared_ptr<kd_tree::KdTree::Node> tree;
struct BinParams {
    float &spheres_dist;
    float &triangles_dist;
    float &cubes_dist;
    Vec3f &hit;
    Vec3f &N;
    raytracing::entities::Material &material;
    BinParams(
    float &spheres_dist,
    float &triangles_dist,
    float &cubes_dist,
    Vec3f &hit,
    Vec3f &N,
    raytracing::entities::Material &material): spheres_dist(spheres_dist), triangles_dist(triangles_dist), cubes_dist(cubes_dist), material(material),
    hit(hit), N(N){}
};

class Render {

public:
    void initialize_kd_tree(std::shared_ptr<raytracing::kd_tree::KdTree::Node> tree);
    static std::vector<raytracing::kd_tree::KdTree::Node *>
    bin_search_in_tree(const raytracing::entities::Ray &ray, std::shared_ptr<raytracing::kd_tree::KdTree::Node> tree);

    /// render - function for filling image with blue color, also for adding spheres and lights to
    ///          picture,
    /// \param spheres - vector of spheres
    /// \param lights - vector of lights
    void render(const char *out_file_path,
                const std::vector<std::unique_ptr<const entities::Figure>> &figures,
                const std::vector<entities::Light> &lights);
};


}//namespace raytracing

#endif //TINYRAYTRACER_RAYTRACING_H
