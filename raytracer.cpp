#include "kd-tree/kdtree.h"
#include "picture/picture.h"
#include "raytracing/raytracing.h"
#include <memory>
#include <variant>

void tree_trace(std::shared_ptr<raytracing::kd_tree::KdTree::Node> tree) {
//    std::cout << "\nOBHOD\n" << std::endl;
    if (tree == nullptr) return;
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


int main(int argc, const char **argv) {
    try
    {
        picture::Picture my_pic(argc, argv);
        raytracing::Render rendering_of_picture;
        rendering_of_picture.initialize_kd_tree(my_pic.kd_tree);
        rendering_of_picture.render(my_pic.out_file_path.c_str(), my_pic.figures, my_pic.lights);
        std::cout << "result is saved in build directory in file " << my_pic.out_file_path.c_str() << std::endl;
    }
    catch (const std::exception &er) {
        std::cout << er.what();
    }
}

