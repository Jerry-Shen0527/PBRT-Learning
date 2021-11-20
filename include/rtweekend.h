
#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <cmath>
#include <limits>
#include <memory>
#include <random>

#ifdef PBRT_Float_AS_Float
typedef Float Float;
#else
typedef float Float;
#endif

// Usings

using std::shared_ptr;
using std::make_shared;
using std::sqrt;

// Constants

const Float infinity = std::numeric_limits<Float>::infinity();
const Float pi = 3.1415926535897932385;

// Utility Functions

inline Float degrees_to_radians(Float degrees) {
    return degrees * pi / 180.0;
}

inline Float random_Float() {
    static std::uniform_real_distribution<Float> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

inline Float random_Float(Float min, Float max) {
    // Returns a random real in [min,max).
    return min + (max - min) * random_Float();
}

inline int random_int(int min, int max) {
    // Returns a random integer in [min,max].
    return static_cast<int>(random_Float(min, max + 1));
}

inline Float clamp(Float x, Float min, Float max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

inline Float Lerp(Float t, Float v1, Float v2) {
    return (1 - t) * v1 + t * v2;
}

template <typename Predicate>
int FindInterval(int size, const Predicate& pred) {
    int first = 0, len = size;
    while (len > 0) {
        int half = len >> 1, middle = first + half;
        // Bisect range based on value of _pred_ at _middle_
        if (pred(middle)) {
            first = middle + 1;
            len -= half + 1;
        }
        else
            len = half;
    }
    return clamp(first - 1, 0, size - 2);
}


#include "ray.h"
#include "vec3.h"

inline vec3 random_cosine_direction() {
    auto r1 = random_Float();
    auto r2 = random_Float();
    auto z = sqrt(1 - r2);

    auto phi = 2 * pi * r1;
    auto x = cos(phi) * sqrt(r2);
    auto y = sin(phi) * sqrt(r2);

    return vec3(x, y, z);
}
#endif