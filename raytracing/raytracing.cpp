#include "raytracing.h"

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <future>

//TODO : description
/// reflect
/// \param I
/// \param N
/// \return
Vec3f reflect(const Vec3f &I, const Vec3f &N) {
    return I - N * 2.f * (I * N);
}

//TODO : description
/// refract
/// \param I
/// \param N
/// \param eta_t
/// \param eta_i
/// \return
Vec3f refract(const Vec3f &I, const Vec3f &N, const float eta_t, const float eta_i=1.f) { // Snell's law
    float cosi = - std::max(-1.f, std::min(1.f, I*N));
    if (cosi<0) return refract(I, -N, eta_i, eta_t); // if the ray comes from the inside the object, swap the air and the media
    float eta = eta_i / eta_t;
    float k = 1 - eta*eta*(1 - cosi*cosi);
    return k<0 ? Vec3f(1,0,0) : I*eta + N*(eta*cosi - sqrtf(k)); // k<0 = total reflection, no ray to refract. I refract it anyways, this has no physical meaning
}

bool PointCheck(const Vec3f A,const Vec3f B,const Vec3f C,const Vec3f CheckPoint)
{
    Vec4f Coord=Vec4f(A.y*(B.z-C.z)+B.y*(C.z-A.z)+C.y*(A.z-B.z),A.z*(B.x-C.x)+B.z*(C.x-A.x)+C.z*(A.x-B.x),A.x*(B.y-C.y)+B.x*(C.y-A.y)+C.x*(A.y-B.y),-(A.x*(B.y*C.z-C.y*B.z)+B.x*(C.y*A.z-A.y*C.z)+C.x*(A.y*B.z-B.y*A.z)));
    return Coord.x * CheckPoint.x + Coord.y * CheckPoint.y + Coord.z * CheckPoint.z + Coord.w == 0;
}


bool scene_intersect(const raytracing::entities::Ray& ray, const std::vector<raytracing::entities::Sphere> &spheres,
                     const std::vector<raytracing::entities::Cube> &cubes, const std::vector<raytracing::entities::Triangle> &triangles,
                     Vec3f &hit, Vec3f &N, raytracing::entities::Material &material) {
    float spheres_dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(ray, dist_i) && dist_i < spheres_dist) {
            spheres_dist = dist_i;
            hit = ray.orig + ray.dir*dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }

    float triangles_dist = std::numeric_limits<float>::max();
    for (size_t i = 0; i < triangles.size(); ++i){
        float dist_i;
        if (triangles[i].ray_intersect(ray, dist_i) && dist_i < triangles_dist){
            triangles_dist = dist_i;
            hit = ray.orig + ray.dir*dist_i;
            N = triangles[i].get_N();
            material = triangles[i].material;
        }
    }

    float cubes_dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < cubes.size(); i++) {
        float dist_i;
        std::vector <float> Coords ={cubes[i].bounds[0].x,cubes[i].bounds[1].x,cubes[i].bounds[0].y,cubes[i].bounds[1].y,cubes[i].bounds[0].z,cubes[i].bounds[1].z};
        float LenD=sqrt(std::pow(Coords[0]-Coords[1],2) + std::pow(Coords[2]-Coords[3],2) + std::pow(Coords[4]-Coords[5],2));
        std::vector <Vec3f> CubePoint = {cubes[i].bounds[0]+Vec3f(0,LenD/sqrt(3),0), cubes[i].bounds[1]+Vec3f(-LenD/sqrt(3),0,0) ,
                                         cubes[i].bounds[1], cubes[i].bounds[1]+Vec3f(0,0,LenD/sqrt(3)), cubes[i].bounds[0],
                                         cubes[i].bounds[0]+Vec3f(0,0,-LenD/sqrt(3)),
                                         cubes[i].bounds[1]+Vec3f(0,-LenD/sqrt(3),0),cubes[i].bounds[0]+Vec3f(LenD/sqrt(3),0,0)};
        if (cubes[i].ray_intersect(ray, dist_i) && dist_i < cubes_dist) {
            cubes_dist = dist_i;
            hit = ray.orig + ray.dir * dist_i;
            Vec3f Vec10 = (CubePoint[0]-CubePoint[1]).normalize();
            Vec3f Vec12 = (CubePoint[2]-CubePoint[1]).normalize();
            Vec3f Vec15 = (CubePoint[5]-CubePoint[1]).normalize();
            Vec3f Vec74 = (CubePoint[4]-CubePoint[7]).normalize();
            Vec3f Vec76 = (CubePoint[6]-CubePoint[7]).normalize();
            Vec3f Vec73 = (CubePoint[3]-CubePoint[7]).normalize();

            if (PointCheck(CubePoint[0],CubePoint[1],CubePoint[2],hit))
                N=cross(Vec10,Vec12);
            if (PointCheck(CubePoint[0],CubePoint[1],CubePoint[5],hit))
                N=cross(Vec10,Vec15);
            if (PointCheck(CubePoint[1],CubePoint[2],CubePoint[5],hit))
                N=cross(Vec12,Vec15);
            if (PointCheck(CubePoint[2],CubePoint[3],CubePoint[6],hit))
                N=cross(Vec73,Vec76);
            if (PointCheck(CubePoint[0],CubePoint[3],CubePoint[4],hit))
                N=cross(Vec73,Vec74);
            if (PointCheck(CubePoint[4],CubePoint[5],CubePoint[6],hit))
                N=cross(Vec74,Vec76);

            //XMinCord=cubes[i].bounds[1].x;

//            std::cout << hit << std::endl;
//            N=Vec3f(0,0,0);
//            if (int(hit.x) == cubes[i].bounds[0].x)
//                N=Vec3f(0,-1,0);
//            if (int(hit.x) == cubes[i].bounds[1].x)
//                N=Vec3f(0,1,0);
//            if (int(hit.y) == cubes[i].bounds[0].y)
//                N=Vec3f(0,0,-1);
//            if (int(hit.y) == cubes[i].bounds[1].y)
//                N=Vec3f(0,0,1);
//            if (int(hit.z) == cubes[i].bounds[0].z)
//                N=Vec3f(-1,0,0);
//            if (int(hit.z) == cubes[i].bounds[1].z)
//                N=Vec3f(1,0,0);
//             TODO вот от N зависит цвет проекции походу каждой
//            N = (hit - cubes[i].bounds[0]).normalize();
            material = cubes[i].material;
        }
    }

    float checkerboard_dist = std::numeric_limits<float>::max();
    if (fabs(ray.dir.y)>1e-3)  {
        float d = -(ray.orig.y+4)/ray.dir.y; // the checkerboard plane has equation y = -4
        Vec3f pt = ray.orig + ray.dir*d;
        if (d>0 && d<spheres_dist && d < triangles_dist) {
            checkerboard_dist = d;
            hit = pt;
            N = Vec3f(0,1,0);
            material.diffuse_color = (int(.5*hit.x+1000) + int(.5*hit.z)) & 1 ? Vec3f(.10, .10, .10) : Vec3f(.3, .2, .1);
        }
    }
    return std::min({spheres_dist, checkerboard_dist, cubes_dist, triangles_dist})<1000;
}


namespace raytracing {

namespace entities{

    Vec3f casting_ray::cast_ray(const Ray &ray,
                                const std::vector<Sphere> &spheres,
                                const std::vector<Light> &lights,
                                const std::vector<entities::Cube> &cubes,
                                const std::vector<entities::Triangle> &triangles, size_t depth) {
        Vec3f point, N;
        Material material;

        if (depth>4 || !scene_intersect(ray, spheres, cubes, triangles, point, N, material)) {
            return Vec3f(0, float(127.0/255), float(255.0/255)); // background color
        }

        ///TODO make just 2 variables
        Vec3f reflect_dir = reflect(ray.dir, N).normalize();
        Vec3f refract_dir = refract(ray.dir, N, material.refractive_index).normalize();
        Ray reflect_ray(reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3,
                                              reflect_dir, &spheres, &lights, &cubes, &triangles, depth+1);// offset the original point to avoid occlusion by the object itself (in first param)
        Ray refract_ray(refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3,
                                              refract_dir, &spheres, &lights, &cubes, &triangles, depth + 1);

        float diffuse_light_intensity = 0, specular_light_intensity = 0;
        for (size_t i=0; i<lights.size(); i++) {
            Vec3f light_dir      = (lights[i].position - point).normalize();
            float light_distance = (lights[i].position - point).norm();

            Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3; // checking if the point lies in the shadow of the lights[i]
            Vec3f shadow_pt, shadow_N;
            Material tmp_material;
            if (scene_intersect(Ray(shadow_orig, light_dir), spheres, cubes, triangles,
                                shadow_pt, shadow_N, tmp_material) &&
                (shadow_pt - shadow_orig).norm() < light_distance)
                continue;

            diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir*N);
            specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*ray.dir),
                                             material.specular_exponent) * lights[i].intensity;
        }
        return material.diffuse_color * diffuse_light_intensity * material.albedo[0] +
               Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1] +
               reflect_ray.color.value()*material.albedo[2] +
               refract_ray.color.value()*material.albedo[3];
    }

    bool Sphere::ray_intersect(const Ray& ray, float &t0) const {
        Vec3f L = center - ray.orig;
        float tca = L*ray.dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }

    //TODO needs testing
    bool Cube::ray_intersect(const Ray &ray, float &t0) const {
        float t_min, t_max, t_y_min, t_y_max, t_z_min, t_z_max;
        if (ray.sign[0]==0 and ray.orig
//        t_min = (bounds[ray.sign[0]].x - ray.orig.x) * ray.invdir.x;
//        t_max = (bounds[1 - ray.sign[0]].x - ray.orig.x) * ray.invdir.x;
//        t_y_min = (bounds[ray.sign[1]].y - ray.orig.y) * ray.invdir.y;
//        t_y_max = (bounds[1 - ray.sign[1]].y - ray.orig.y) * ray.invdir.y;
//        t_z_min = (bounds[ray.sign[2]].z - ray.orig.z) * ray.invdir.z;
//        t_z_max = (bounds[1-ray.sign[2]].z - ray.orig.z) * ray.invdir.z;
//        if ((t_min > t_y_max) || (t_y_min > t_max))
//            return false;
//        if ((t_y_min > t_z_max) || (t_z_min > t_y_max))
//            return false;
//        if (t_y_min > t_min)
//            t_min = t_y_min;
//        if (t_y_max < t_max)
//            t_max = t_y_max;
//
//        t_z_min = (bounds[ray.sign[2]].z - ray.orig.z) * ray.invdir.z;
//        t_z_max = (bounds[1 - ray.sign[2]].z - ray.orig.z) * ray.invdir.z;
//
//        if ((t_min > t_z_max) || (t_z_min > t_max))
//            return false;
//        if (t_z_min > t_min)
//            t_min = t_z_min;
//        if (t_z_max < t_max)
//            t_max = t_z_max;
//
//        t0 = t_min;
//
//        if (t0 < 0) {
//            t0 = t_max;
//            if (t0 < 0) return false;
//        }
//
//        return true;
    }

    bool Triangle::ray_intersect(const Ray &ray, float &t0) const {
        const float EPSILON = 0.0000001;
        Vec3f vertex0 = p0;
        Vec3f vertex1 = p1;
        Vec3f vertex2 = p2;
        Vec3f edge1, edge2, h, s, q;
        float a, f, u, v;
        edge1 = vertex1 - vertex0;
        edge2 = vertex2 - vertex0;
        h = cross(ray.dir, edge2);
        a = edge1*h;
        if (a > -EPSILON && a < EPSILON)
            return false; // This ray is parallel to this triangle.
        f = 1.0 / a;
        s = ray.orig - vertex0;
        u = f * (s * h);
        if (u < 0.0 || u > 1.0)
            return false;
        q = cross(s, edge1);
        v = f * (ray.dir * q);
        if (v < 0.0 || u + v > 1.0)
            return false;
        // At this stage we can compute t to find out where the intersection point is on the line.
        float t = f * (edge2 * q);
        if (t > EPSILON) // ray intersection
        {
            t0 = t;
            return true;
        }
        else { // This means that there is a line intersection but not a ray intersection.
            return false;
        }
    }
}// namespace entities




void render(const char* out_file_path, const std::vector<entities::Sphere> &spheres, const std::vector<entities::Light> &lights,
            const std::vector<entities::Cube> &cubes, const std::vector<entities::Triangle> &triangles) {
    const int   width    = 1920;
    const int   height   = 1080;
    const float fov      = M_PI/3.; ///that's a viewing angle = pi/3
    std::vector<Vec3f> framebuffer(width*height);
    const auto amount_of_threads = 1;//std::thread::hardware_concurrency();
    std::vector<std::future<void>> tasks(amount_of_threads);

    size_t portion = width / amount_of_threads;
    for (size_t j = 0; j<height; j++) { // actual rendering loop
        //trying to parallel compute pieces of "line" of a picture
        for (size_t start = 0, finish = portion, index = 0;
        index < amount_of_threads; ++index, start=finish,
        finish = index == amount_of_threads - 1 ? width : finish + portion){
            tasks[index] = std::async(std::launch::async, [=, &spheres, &lights, &cubes, &triangles, &framebuffer](){
                for (size_t i = start; i < finish; i++) {
                    auto dir_x = (i + 0.5) - width / 2.;
                    auto dir_y = -(j + 0.5) + height / 2.;    // this flips the image at the same time
                    auto dir_z = -height / (2. * tan(fov / 2.));
                    framebuffer[i + j * width] = entities::casting_ray::cast_ray(entities::Ray(
                            Vec3f(0, 0, 0), Vec3f(dir_x, dir_y, dir_z).normalize()),
                                                          spheres, lights, cubes, triangles);
                }
            });
        }
        /// waiting our parallel tasks
        for (auto&& task : tasks){
            task.get();
        }
    }

    std::ofstream ofs; // save the framebuffer to file
    ofs.open(out_file_path, std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c*(1./max);
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

}// namespace raytracing
