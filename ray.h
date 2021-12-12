#ifndef RAY_H
#define RAY_H

#include "vec3.h"
//struct hit_triangle
//{
//    bool hitted_triangle = false;
//    int kx, ky, kz;
//    vec3 d;
//    double Sx, Sy, Sz;
//
//
//    inline void setxyz(const ray& r)
//    {
//        hitted_triangle = true;
//        kz = maxdimension(vector_abs(r.dir));
//        kx = kz + 1;
//        if (kx == 3) kx = 0;
//        ky = kx + 1;
//        if (ky == 3)ky = 0;
//        d = permutes(r.dir, kx, ky, kz);
//        double Sx = -d.e[0] / d.e[2];
//        double Sy = -d.e[1] / d.e[2];
//        double Sz = 1.f / d.e[2];
//    }
//};

class ray {
public:
    ray() {}
    ray(const point3& origin, const vec3& direction, double time = 0.0)
        : orig(origin), dir(direction), tm(time)
    {}

    point3 origin() const { return orig; }
    vec3 direction() const { return dir; }
    double time() const { return tm; }
    //void setxyz();
    point3 at(double t) const {
        return orig + t * dir;
    }

public:
    point3 orig;
    vec3 dir;
    double tm;


    bool hitted_triangle = false;
    int kx, ky, kz;
    vec3 d;
    double Sx, Sy, Sz;


};

//void ray::setxyz()
//{
//    hitted_triangle = true;
//    kz = maxdimension(vector_abs(r.dir));
//    kx = kz + 1;
//    if (kx == 3) kx = 0;
//    ky = kx + 1;
//    if (ky == 3)ky = 0;
//    d = permutes(r.dir, kx, ky, kz);
//    double Sx = -d.e[0] / d.e[2];
//    double Sy = -d.e[1] / d.e[2];
//    double Sz = 1.f / d.e[2];
//}
#endif