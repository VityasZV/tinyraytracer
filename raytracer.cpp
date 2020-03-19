#include "picture/picture.h"
#include "raytracing/raytracing.h"

int picture::Picture::scene_id=2;
int main(int argc, const char **argv) {
    try
    {
        picture::Picture my_pic(argc, argv);
        raytracing::render(my_pic.out_file_path.c_str(), my_pic.figures, my_pic.lights);
        std::cout << "result is saved in build directory in file " << my_pic.out_file_path.c_str() << std::endl;
    }
    catch (const std::exception &er) {
        std::cout << er.what();
    }
    return 0;
}

