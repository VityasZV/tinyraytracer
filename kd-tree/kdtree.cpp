//
// Created by Виктор Задябин on 16.03.2020.
//

#include "kdtree.h"

#include <utility>
#include <future>
#include <tuple>

namespace raytracing::kd_tree {

std::shared_ptr<KdTree::Node>
KdTree::build(const std::vector<std::shared_ptr<RenderWrapper>> &objs, entities::AABB &box, int depth) {
    if (objs.empty() or depth > KdTree::max_depth) {
        return nullptr;
    }
    if (objs.size() <= KdTree::min_amount_of_figures) {
        std::variant<raytracing::entities::Plane, const std::vector<std::shared_ptr<RenderWrapper>>> var = objs;
        return std::make_shared<KdTree::Node>(box, objs);
    }
    //TODO make SAH algorithm because thats shit
    // algorithm 1 (naive kdtree)
    auto best_plane = cut(objs, entities::Axis::x); //fixed plane on x axis
    std::optional<std::pair<raytracing::entities::AABB, raytracing::entities::AABB>> pair_of_cubes;
    try {
        pair_of_cubes = box.cut(best_plane);
    } catch (const std::exception &er) {
        //добавить объекты
        std::variant<raytracing::entities::Plane, const std::vector<std::shared_ptr<RenderWrapper>>> var = objs;
        //plane is outside of the box
        return std::make_shared<KdTree::Node>(box, var);
    }
    //иначе добавим плоскость
    std::vector<std::shared_ptr<RenderWrapper>> objl, objr;
    for (auto &obj : objs) {
        if (obj->box.GetVMax()[best_plane.GetAxis()] >= best_plane.GetPos() - EPS) {
            objr.push_back(obj);
        }
        if (obj->box.GetVMin()[best_plane.GetAxis()] <= best_plane.GetPos() + EPS) {
            objl.push_back(obj);
        }
    }
    std::variant<raytracing::entities::Plane, const std::vector<std::shared_ptr<RenderWrapper>>> var = best_plane;

    if (depth < std::thread::hardware_concurrency()) {
        auto left_node_ptr = std::async(std::launch::async, [&]() {
            return build(objl, pair_of_cubes->first, depth + 1);
        });
        auto right_node_ptr = std::async(std::launch::async, [&]() {
            return build(objr, pair_of_cubes->second, depth + 1);
        });
        return std::make_shared<KdTree::Node>(box, var, left_node_ptr.get(), right_node_ptr.get());
    } else {
        return std::make_shared<KdTree::Node>(box, var, build(objl, pair_of_cubes->first, depth + 1),
                                              build(objr, pair_of_cubes->second, depth + 1));
    }

}

raytracing::entities::Plane
KdTree::cut(const std::vector<std::shared_ptr<RenderWrapper>> &objs, const raytracing::entities::Axis &axis) {
    std::vector<float> min_list;
    min_list.resize(objs.size());
    size_t i = 0;
    for (auto &obj_ptr : objs) {
        min_list[i++] = obj_ptr->box.GetVMin()[axis];
    }
    //sort that relates to [min_list.begin() + min_list.size() / 2] element
    std::nth_element(min_list.begin(), min_list.begin() + min_list.size() / 2, min_list.end());
    raytracing::entities::Plane ret(axis, min_list[min_list.size() / 2] + 2 * EPS);
    return ret;
}

//std::shared_ptr<Trace> KdTree::get_trace(const raytracing::entities::Ray &ray, float max_dist) const {
//
//}


}// namespace raytracing::kd_tree