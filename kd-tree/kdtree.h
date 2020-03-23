//
// Created by Виктор Задябин on 16.03.2020.
//

#ifndef RAYTRACER_KDTREE_H
#define RAYTRACER_KDTREE_H

#include <variant>
#include "../raytracing/entities.h"

namespace picture {
class Picture;
}
namespace raytracing::kd_tree {
enum class NodeOrLeaf {
    Node, Leaf
};

class KdTree {
public:
    const static unsigned int min_amount_of_figures = 2;
    const static unsigned int max_depth = 20;
    struct RenderWrapper;

    struct Node {
        raytracing::entities::AABB box;
        std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>> child;
        //NodeOrLeaf node_or_leaf;
        //Plane for non-leaf node and figures for leaf node
        std::variant<raytracing::entities::Plane, const std::vector<std::shared_ptr<RenderWrapper>>> plane_or_figures;

        Node(const raytracing::entities::AABB &box,
             std::variant<raytracing::entities::Plane, const std::vector<std::shared_ptr<RenderWrapper>>> plane_or_figures,
             std::shared_ptr<Node> p1 = nullptr,
             std::shared_ptr<Node> p2 = nullptr)
                : box(box), plane_or_figures(std::move(plane_or_figures)) {
            if (p1 == nullptr && p2 == nullptr) {
                //node_or_leaf = NodeOrLeaf::Leaf;

            } else {
                // node_or_leaf = NodeOrLeaf::Node;
                child.first = std::move(p1);
                child.second = std::move(p2);
            }
        }

        ~Node() = default;
    };

    struct RenderWrapper {
        std::shared_ptr<const raytracing::entities::Figure> obj;
        raytracing::entities::AABB box;

        RenderWrapper(std::shared_ptr<const raytracing::entities::Figure> obj_p, const raytracing::entities::AABB &box)
                : box(box), obj(obj_p) {}
    };

    //std::shared_ptr<Trace> get_trace(const raytracing::entities::Ray &ray, float max_dist) const;

private:
    static std::shared_ptr<Node>
    build(const std::vector<std::shared_ptr<RenderWrapper>> &objs, raytracing::entities::AABB &box, int depth);

    [[nodiscard]] static raytracing::entities::Plane
    cut(const std::vector<std::shared_ptr<RenderWrapper>> &objs, const raytracing::entities::Axis &axis);

    friend class picture::Picture;
};

}// namespace raytracing::kd_tree


#endif //RAYTRACER_KDTREE_H
