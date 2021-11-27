
#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>
#include "rtweekend.h"

//#define PBRT_DEBUG
#ifndef PBRT_DEBUG

#define IN_RANGE(condition) (void)0
#define ZERO_DENOMINATOR(t) (void)0

#else
#define IN_RANGE(condition) if (!condition) std::cerr << "vec:Out of range" << std::endl;
#define ZERO_DENOMINATOR(t) if (abs(t)<1e-8) std::cerr << "denominator is near zero." << std::endl;

#endif // !PBRT_DEBUG


using std::sqrt;

class vec3 {
public:
    vec3() { x = y = z = 0; }
    vec3(Float e0, Float e1, Float e2) : x(e0),y(e1),z(e2) {}

    

    vec3 operator-() const { return vec3(-x, -y, -z); }
    Float operator[](int i) const {
        IN_RANGE((i >= 0 && i < 3));
        if (i == 0) return x;
        if (i == 1) return y;
        if (i == 2) return z;

    }
    Float& operator[](int i) { 
        IN_RANGE((i >= 0 && i < 3));
        if (i == 0) return x;
        if (i == 1) return y;
        if (i == 2) return z;
    }

    vec3& operator+=(const vec3& v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    vec3& operator-=(const vec3& v)
    {
        return *this += (-v);
    }

    vec3& operator*=(const Float t) {
        x *= t;
        y *= t;
        z *= t;
        return *this;
    }

    vec3& operator/=(const Float t) {
        ZERO_DENOMINATOR(t);
        return *this *= 1 / t;
    }

    Float length() const {
        return sqrt(length_squared());
    }

    Float length_squared() const {
        return x * x + y * y + z * z;
    }

    inline static vec3 random() {
        return vec3(random_Float(), random_Float(), random_Float());
    }

    inline static vec3 random(Float min, Float max) {
        return vec3(random_Float(min, max), random_Float(min, max), random_Float(min, max));
    }

    bool near_zero() const {
        // Return true if the vector is close to zero in all dimensions.
        const auto s = 1e-8;
        return (fabs(x) < s) && (fabs(y) < s) && (fabs(z) < s);
    }
   /*Float x() const { return x; }
    Float y() const { return y; }
    Float z() const { return z; }*/
public:
    Float x, y, z;
};


// Type aliases for vec3
using point3 = vec3;   // 3D point
using color = vec3;    // RGB color
using Normal = vec3;
// vec3 Utility Functions

inline std::ostream& operator<<(std::ostream& out, const vec3& v) {
    return out << v.x << ' ' << v.y << ' ' << v.z;
}

inline vec3 operator+(const vec3& u, const vec3& v) {
    return vec3(u.x + v.x, u.y + v.y, u.z + v.z);
}

inline vec3 operator-(const vec3& u, const vec3& v) {
    return vec3(u.x - v.x, u.y - v.y, u.z - v.z);
}

inline vec3 operator*(const vec3& u, const vec3& v) {
    return vec3(u.x * v.x, u.y * v.y, u.z * v.z);
}

inline vec3 operator*(Float t, const vec3& v) {
    return vec3(t * v.x, t * v.y, t * v.z);
}

inline vec3 operator*(const vec3& v, Float t) {
    return t * v;
}

inline vec3 operator/(vec3 v, Float t) {
    return (1 / t) * v;
}

inline Float dot(const vec3& u, const vec3& v) {
    return u.x * v.x
        + u.y * v.y
        + u.z * v.z;
}

inline vec3 cross(const vec3& u, const vec3& v) {
    return vec3(u.y * v.z - u.z * v.y,
        u.z * v.x - u.x * v.z,
        u.x * v.y - u.y * v.x);
}

inline Float Distance(const vec3& u, const vec3& v)
{
    return (u - v).length();
}

inline vec3 unit_vector(vec3 v) {
    return v / v.length();
}

inline Normal Faceforward(const Normal& n, const vec3& v)
{
    return dot(n, v) >= 0.f ? n : -n;
}

inline vec3 random_in_unit_sphere() {
        while (true) {
            auto p = vec3::random(-1, 1);
            if (p.length_squared() >= 1) continue;
            return p;
        }
    }

inline vec3 random_unit_vector() {
     return unit_vector(random_in_unit_sphere());
 }

 inline vec3 random_in_hemisphere(const vec3& normal) {
     vec3 in_unit_sphere = random_in_unit_sphere();
     if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
         return in_unit_sphere;
     else
         return -in_unit_sphere;
 }
 
 inline vec3 reflect(const vec3& v, const vec3& n) {
     return v - 2 * dot(v, n) * n;
 }

 inline vec3 refract(const vec3& uv, const vec3& n, Float etai_over_etat) {
     auto cos_theta = fmin(dot(-uv, n), 1.0);
     vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
     vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
     return r_out_perp + r_out_parallel;
 }

 inline vec3 random_in_unit_disk() {
     while (true) {
         auto p = vec3(random_Float(-1, 1), random_Float(-1, 1), 0);
         if (p.length_squared() >= 1) continue;
         return p;
     }
 }
#endif