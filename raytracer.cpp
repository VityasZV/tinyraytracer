#include "kd-tree/kdtree.h"
#include "picture/picture.h"
#include "raytracing/raytracing.h"
#include <memory>
#include <variant>



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

