#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "picture/picture.h"
#include "raytracing/raytracing.h"


int main(int argc, const char **argv) {
    try {
        int n = -1;
        int envmap_width, env_height;
        unsigned char *pixmap = stbi_load("envmap.jpg", &envmap_width, &env_height, &n, 0);
        picture::Picture my_pic(argc, argv);
        switch (my_pic.scene_id) {
            case 1:
                raytracing::render(my_pic.out_file_path.c_str(), my_pic.figures, my_pic.lights);
                std::cout << "result is saved in build directory in file " << my_pic.out_file_path.c_str() << std::endl;
                break;
            case 2:
                return 0;
            case 3:
                return 0;
            default:
                throw std::runtime_error("Incorrect scene number!");
        }

    }
    catch (const std::exception &er) {
        std::cout << er.what();
    }
    return 0;
}

