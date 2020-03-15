//
// Created by Виктор Задябин on 14.03.2020.
//

#include "picture.h"

#include <fstream>
#include <sstream>
#include <future>

picture::Picture::Picture(const int argc, const char **argv) {
    PreparingOutFileAndScene(argc, argv);
    for (const auto &p : spheres_params) {
        figures.emplace_back(std::make_unique<raytracing::entities::Sphere>(p.coordinates, p.radius, p.material));
    }
    for (const auto &p : lights_params) {
        lights.emplace_back(raytracing::entities::Light(p.position, p.intensity));
    }
    const auto shift1 = Vec3f(5, 1, -6);
    const auto shift2 = Vec3f(-5, 1, -6);
    const auto shift3 = Vec3f(15, 1, -6);
    const auto shifts = {shift1, shift2, shift3};
    //here is my triangles
    for (const auto &p : triangle_params) {
        for (const auto &s : shifts) {
            figures.emplace_back(
                    std::make_unique<raytracing::entities::Triangle>(p.p0 + s, p.p1 + s, p.p2 + s,
                                                                     p.material));
        }
    }
    triangle_params.clear();
    //here comes the duck
    MakeTriangleMash("../duck.obj");
    //here comes the deer
    MakeTriangleMash("../deer.obj");
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
}

void picture::Picture::MakeTriangleMash(const char *file_name) {
    auto shift_deer = Vec3f(-5, -4, -8);
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
                verticels.push_back(v);
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


