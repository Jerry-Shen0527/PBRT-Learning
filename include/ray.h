#ifndef RAY_H
#define RAY_H

#include "vec3.h"

class ray {
public:
    ray():mint(0.f),maxt(infinity),depth(0){}
    ray(const point3& origin, const vec3& direction, Float time = 0.0,Float min_t=0.f,Float max_t=infinity,int depth=0)
        : orig(origin), dir(direction), tm(time), mint(min_t), maxt(max_t), depth(depth)
    {}
    ray(const point3& origin, const vec3& direction, const ray& parent, Float min_t = 0.f, Float max_t = infinity)
        :orig(origin),dir(direction),tm(parent.tm), mint(min_t), maxt(max_t),depth(parent.depth+1)
    {}

    point3 origin() const { return orig; }
    vec3 direction() const { return dir; }
    Float time() const { return tm; }

    point3 at(Float t) const {
        return orig + t * dir;
    }

public:
    point3 orig;
    vec3 dir;
    Float tm;
    mutable Float mint, maxt;
    int depth;
};

#endif
