
#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <cmath>
#include <limits>
#include <memory>
#include <random>

#ifdef PBRT_Float_AS_Float
typedef double Float;
#else
typedef float Float;
#endif

#ifdef _MSC_VER
#define MachineEpsilon (std::numeric_limits<Float>::epsilon() * 0.5)
#else
static PBRT_CONSTEXPR Float MachineEpsilon =
std::numeric_limits<Float>::epsilon() * 0.5;
#endif

//#define PBRT_DEBUG
#ifndef PBRT_DEBUG

#define IN_RANGE(condition) (void)0
#define ZERO_DENOMINATOR(t) (void)0
#define CHECK_LE(t1,t2) (void)0
#define CHECK_LT(t1,t2) (void)0
#define CHECK_GE(t1,t2) (void)0
#define CHECK_GT(t1,t2) (void)0
#define CHECK_NE(t1,t2) (void)0
#define CHECK_EQ(t1,t2) (void)0

#else
#include <iostream>
#define IN_RANGE(condition) if (!condition) std::cerr << "vec:Out of range" << std::endl;
#define ZERO_DENOMINATOR(t) if (t<1e-8 && t>-1e-8) std::cerr << "denominator is near zero." << std::endl;
#define CHECK_LE(t1,t2) if(t1 > t2) std::cerr<< "Not less/equal than" << std::endl;
#define CHECK_LT(t1,t2) if(t1 >= t2) std::cerr<< "Not less than" << std::endl;
#define CHECK_GE(t1,t2) if(t1 < t2) std::cerr<< "Not greater/equal than" <<std::endl;
#define CHECK_GT(t1,t2) if(t1 <= t2) std::cerr<< "Not greater than" <<std::endl;
#define CHECK_NE(t1,t2) if(t1 == t2) std::cerr<<t1<<" equals to "<<t2<<std::endl;
#define CHECK_EQ(t1,t2) if(t1 != t2) std::cerr<<t1<<" not equals to "<<t2<<std::endl;

#endif // !PBRT_DEBUG


// Usings

using std::shared_ptr;
using std::make_shared;
using std::sqrt;

// Constants

const Float Infinity = std::numeric_limits<Float>::infinity();
const Float Pi = 3.1415926535897932385;

// Utility Functions
inline uint32_t FloatToBits(float f) {
    uint32_t ui;
    memcpy(&ui, &f, sizeof(float));
    return ui;
}

inline float BitsToFloat(uint32_t ui) {
    float f;
    memcpy(&f, &ui, sizeof(uint32_t));
    return f;
}

inline uint64_t FloatToBits(double f) {
    uint64_t ui;
    memcpy(&ui, &f, sizeof(double));
    return ui;
}

inline double BitsToFloat(uint64_t ui) {
    double f;
    memcpy(&f, &ui, sizeof(uint64_t));
    return f;
}

inline float NextFloatUp(float v) {
    // Handle infinity and negative zero for _NextFloatUp()_
    if (std::isinf(v) && v > 0.) return v;
    if (v == -0.f) v = 0.f;

    // Advance _v_ to next higher float
    uint32_t ui = FloatToBits(v);
    if (v >= 0)
        ++ui;
    else
        --ui;
    return BitsToFloat(ui);
}

inline float NextFloatDown(float v) {
    // Handle infinity and positive zero for _NextFloatDown()_
    if (std::isinf(v) && v < 0.) return v;
    if (v == 0.f) v = -0.f;
    uint32_t ui = FloatToBits(v);
    if (v > 0)
        --ui;
    else
        ++ui;
    return BitsToFloat(ui);
}

inline double NextFloatUp(double v, int delta = 1) {
    if (std::isinf(v) && v > 0.) return v;
    if (v == -0.f) v = 0.f;
    uint64_t ui = FloatToBits(v);
    if (v >= 0.)
        ui += delta;
    else
        ui -= delta;
    return BitsToFloat(ui);
}

inline double NextFloatDown(double v, int delta = 1) {
    if (std::isinf(v) && v < 0.) return v;
    if (v == 0.f) v = -0.f;
    uint64_t ui = FloatToBits(v);
    if (v > 0.)
        ui -= delta;
    else
        ui += delta;
    return BitsToFloat(ui);
}

inline Float gamma(int n) {
    return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}

inline Float Abs(Float x) {
    return x >= 0 ? x : -x;
}

inline Float Radians(Float deg) { return (Pi / Float(180)) * deg; }

inline Float Degrees(Float rad) { return (Float(180) / Pi) * rad; }


inline Float RandomFloat() {
    static std::uniform_real_distribution<Float> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

inline Float RandomFloat(Float min, Float max) {
    // Returns a random real in [min,max).
    return min + (max - min) * RandomFloat();
}

inline int RandomInt(int min, int max) {
    // Returns a random integer in [min,max].
    return static_cast<int>(RandomFloat(min, max + 1));
}

inline Float Clamp(Float x, Float min, Float max) {
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
    return Clamp(first - 1, 0, size - 2);
}





#endif