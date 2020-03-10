#include "picture.h"
#include "raytracing/raytracing.h"

int main() {
    try {
        picture::Picture my_pic;
        raytracing::render(my_pic.spheres, my_pic.lights);
        std::cout << "result is saved in build directory in file out.ppm" << std::endl;
    }
    catch (const std::exception& er){
        std::cout << er.what();
    }
    return 0;
}

