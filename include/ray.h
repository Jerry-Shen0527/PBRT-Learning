#ifndef RAY_H
#define RAY_H

#include "vec3.h"

class ray {
public:
    ray() {}
    ray(const point3& origin, const vec3& direction, Float time = 0.0)
        : orig(origin), dir(direction), tm(time)
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
};

#endif
