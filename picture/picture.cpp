//
// Created by Виктор Задябин on 14.03.2020.
//

#include "picture.h"

#include <fstream>
#include <sstream>
#include <future>
#include <cstring>
#include <string.h>


picture::Picture::  Picture(const int argc, const char **argv) {
    PreparingOutFileAndScene(argc, argv);
//    for (const auto &p : spheres_params) {
//        figures.emplace_back(std::make_unique<raytracing::entities::Sphere>(p.coordinates, p.radius, p.material));
//    }
    for (const auto &p : lights_params) {
        lights.emplace_back(raytracing::entities::Light(p.position, p.intensity));
    }
    const auto shift1 = Vec3f(5, 1, -6);
    const auto shift2 = Vec3f(-5, 1, -6);
    const auto shift3 = Vec3f(15, 1, -6);
    const auto shifts = {shift1, shift2, shift3};
    //here is my triangles
//    for (const auto &p : triangle_params) {
//        for (const auto &s : shifts) {
//            figures.emplace_back(
//                    std::make_unique<raytracing::entities::Triangle>(p.p0 + s, p.p1 + s, p.p2 + s,
//                                                                     p.material));
//        }
//    }
    triangle_params.clear();
    //here comes the duck
    const auto duck_shift1 = Vec3f{-10, 0, 0};
    const auto duck_shift2 = Vec3f{-10, 0, -10};
    const auto duck_shift3 = Vec3f{0, 0, -10};


    MakeTriangleMash("../duck.obj", Vec3f(0,0,0));
    MakeTriangleMash("../duck.obj", duck_shift1);
    MakeTriangleMash("../duck.obj", duck_shift2);
    MakeTriangleMash("../duck.obj", duck_shift3);
    //here comes the deer
    //MakeTriangleMash("../deer.obj");
    //std::cout << "Всего примитивов " << figures.size() << std::endl;
    FormKdTree();
}

void picture::Picture::PreparingOutFileAndScene(int argc, const char **argv) {
    std::unordered_map<std::string, std::string> cmd_line_params;
    for (int i = 0; i < argc; ++i) {
        std::string key(argv[i]);
        if (key.size() > 0 && key[0] == '-') {
            if (i != argc - 1) // not last argument
            {
                cmd_line_params[key] = argv[i + 1];
                ++i;
            } else
                cmd_line_params[key] = "";
        }
    }
    if (cmd_line_params.find("-out") != cmd_line_params.end())
        out_file_path = cmd_line_params["-out"];
    if (cmd_line_params.find("-scene") != cmd_line_params.end())
        scene_id = atoi(cmd_line_params["-scene"].c_str());
    switch (scene_id) {
        case 1:
            spheres_params = {
                    {Vec3f(-3, 0, -16),      2, Materials[MaterialName::ivory]},
                    {Vec3f(-1.0, -1.5, -12), 2, Materials[MaterialName::glass]},
                    {Vec3f(1.5, -0.5, -18),  3, Materials[MaterialName::red_rubber]},
                    {Vec3f(7, 5, -18),       4, Materials[MaterialName::mirror]},
                    {Vec3f(-8, 5, -18),      4, Materials[MaterialName::mirror]}
            };
            triangle_params = {
                    //dont forget about right trio while adding params!!!
                    //front
                    {Vec3f(-6, 0, -6), Vec3f(-5, 0, -6), Vec3f(-5, 2, -6), Materials[MaterialName::red_rubber]},
                    {Vec3f(-5, 2, -6), Vec3f(-6, 2, -6), Vec3f(-6, 0, -6), Materials[MaterialName::red_rubber]},
                    //back
                    {Vec3f(-5, 0, -6), Vec3f(-6, 0, -6), Vec3f(-5, 0, -8), Materials[MaterialName::red_rubber]},
                    {Vec3f(-5, 0, -8), Vec3f(-6, 0, -6), Vec3f(-6, 0, -8), Materials[MaterialName::red_rubber]},
                    //up
                    {Vec3f(-6, 2, -6), Vec3f(-5, 2, -6), Vec3f(-5, 2, -8), Materials[MaterialName::red_rubber]},
                    {Vec3f(-5, 2, -8), Vec3f(-6, 2, -8), Vec3f(-6, 2, -6), Materials[MaterialName::red_rubber]},
                    //left
                    {Vec3f(-6, 2, -8), Vec3f(-6, 2, -6), Vec3f(-6, 0, -6), Materials[MaterialName::red_rubber]},
                    {Vec3f(-6, 0, -8), Vec3f(-6, 2, -8), Vec3f(-6, 0, -6), Materials[MaterialName::red_rubber]},
                    //right
                    {Vec3f(-5, 0, -8), Vec3f(-5, 2, -8), Vec3f(-5, 2, -6), Materials[MaterialName::red_rubber]},
                    {Vec3f(-5, 0, -8), Vec3f(-5, 2, -6), Vec3f(-5, 0, -6), Materials[MaterialName::red_rubber]},
                    //down
                    {Vec3f(-6, 2, -8), Vec3f(-6, 0, -8), Vec3f(-5, 0, -8), Materials[MaterialName::red_rubber]},
                    {Vec3f(-6, 2, -8), Vec3f(-5, 2, -8), Vec3f(-5, 0, -8), Materials[MaterialName::red_rubber]},

            };
            break;
        case 2:
            spheres_params = {
                    {Vec3f(-10, -2.5, -10),    2, Materials[MaterialName::red_rubber]},
                    {Vec3f(-10, -2.5, -30), 2, Materials[MaterialName::red_rubber]},
                    {Vec3f(10, -2.5, -10),  2, Materials[MaterialName::red_rubber]},
                    {Vec3f(10, -2.5, -30),      2, Materials[MaterialName::red_rubber]},
                    {Vec3f(-8, 5, -18),     4, Materials[MaterialName::mirror]},
                    {Vec3f(7, 5, -18),      4, Materials[MaterialName::mirror]}

            };
            break;
        case 3:
        default:
            throw std::runtime_error("Incorrect scene number!");
    }
}

void picture::Picture::MakeTriangleMash(const char *file_name, const Vec3f& shift) {
    auto shift_deer = Vec3f(-5, -4, -8);
    shift_deer = shift_deer + shift;
    auto figure_material = Materials[MaterialName::red_rubber];
    std::vector<Vec3f> verticels;
    std::vector<Vec2f> uvIndices;
    std::vector<Vec3i> faces;
    std::ifstream in(file_name, std::ifstream::in);
    if (in.fail()) {
        std::stringstream error;
        error << "Failed to open " << file_name;
        throw std::runtime_error(error.str());
    }
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i = 0; i < 3; ++i) iss >> v[i];
            if (strncmp(file_name, "../deer.obj", 11) == 0) {
                verticels.push_back(v / 200. + shift_deer);
            } else {
                verticels.push_back(v + shift);
            }
        } else if (!line.compare(0, 2, "vt")) {
            continue;
        } else if (!line.compare(0, 2, "f ")) {
            if (std::string(file_name) == std::string("../duck.obj")) {
                Vec3i f;
                int idx, cnt = 0;
                iss >> trash;
                while (iss >> idx) {
                    --idx; // in wavefront obj all indices start at 1, not zero
                    f[cnt++] = idx;
                }
                if (3 == cnt) faces.push_back(f);
            } else {
                Vec3i f;
                int idx, cnt = 0;
                iss >> trash;
                while (iss >> idx) {
                    --idx; // in wavefront obj all indices start at 1, not zero
                    f[cnt++] = idx;
                    while (iss.get() != ' ') {}
                }
                if (3 == cnt) faces.push_back(f);
            }
        }
    }
    for (const auto &face : faces) {
        figures.emplace_back(std::make_unique<raytracing::entities::Triangle>(verticels[face.x], verticels[face.y],
                                                                              verticels[face.z], figure_material));
    }

}

void picture::Picture::FormKdTree() {
    raytracing::entities::AABB space(Vec3f(-30, -16, -5), Vec3f(30, 16, -30)); //initial coube
    //DONE
    //пихать в figures_in_tree_root нужно не внешний куб, а куб для каждого блять примитива на хуй -- done
    std::vector<std::shared_ptr<raytracing::kd_tree::KdTree::RenderWrapper>> figures_in_tree_root;
    for (auto &figure : figures) {
        auto bounding_box = figure->GetAABB();
        figures_in_tree_root.push_back(
                std::make_shared<raytracing::kd_tree::KdTree::RenderWrapper>(std::move(figure), bounding_box));
    }
    figures.clear();
    kd_tree = raytracing::kd_tree::KdTree::build(figures_in_tree_root, space, 0);
    //std::cout << "finish" << std::endl;
}


