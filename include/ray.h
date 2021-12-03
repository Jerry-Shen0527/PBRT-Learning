#ifndef RAY_H
#define RAY_H

#include "vec3.h"

class ray {
public:
    ray():tMin(0.f),tMax(Infinity),depth(0){}
    ray(const point3& origin, const vec3& direction, Float time = 0.0,Float min_t=0.f,Float max_t=Infinity,int depth=0)
        : o(origin), d(direction), time(time), tMin(min_t), tMax(max_t), depth(depth)
    {}
    ray(const point3& origin, const vec3& direction, const ray& parent, Float min_t = 0.f, Float max_t = Infinity)
        :o(origin),d(direction),time(parent.time), tMin(min_t), tMax(max_t),depth(parent.depth+1)
    {}
    Point3f operator()(Float t) const { return o + d * t; }

    point3 origin() const { return o; }
    vec3 direction() const { return d; }
    Float Time() const { return time; }

    point3 at(Float t) const {
        return o + t * d;
    }

public:
    point3 o;
    vec3 d;
    Float time;
    mutable Float tMin, tMax;
    int depth;
};


using Ray = ray;

#endif
