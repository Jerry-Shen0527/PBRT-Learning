#pragma once
#include <cmath>
#include <iostream>
#include "random.h"
const double pi = 3.1415926535897932385;
using std::sqrt;

class Vector3f {
public:
    Vector3f() : e{ 0,0,0 } {}
    Vector3f(double e0, double e1, double e2) : e{ e0, e1, e2 } {}

    double x() const { return e[0]; }
    double y() const { return e[1]; }
    double z() const { return e[2]; }

    Vector3f operator-() const { return Vector3f(-e[0], -e[1], -e[2]); }
    double operator[](int i) const { return e[i]; }
    double& operator[](int i) { return e[i]; }

    Vector3f& operator+=(const Vector3f& v) {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    Vector3f& operator*=(const double t) {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    Vector3f& operator/=(const double t) {
        return *this *= 1 / t;
    }

    double length() const {
        return sqrt(length_squared());
    }

    double length_squared() const {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }

    inline static Vector3f random() {
        return Vector3f(rand() / (RAND_MAX + 1.0), rand() / (RAND_MAX + 1.0), rand() / (RAND_MAX + 1.0));
    }

    inline static Vector3f random(double min, double max) {
        return Vector3f(min + (max - min) * rand() / (RAND_MAX + 1.0),
            min + (max - min) * rand() / (RAND_MAX + 1.0),
            min + (max - min) * rand() / (RAND_MAX + 1.0));
    }

    bool near_zero() const {
        // Return true if the vector is close to zero in all dimensions.
        const auto s = 1e-8;
        return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
    }


public:
    double e[3];
};

// Type aliases for Vector3f
using Point3f = Vector3f;   // 3D point
using color = Vector3f;    // RGB color

inline std::ostream& operator<<(std::ostream& out, const Vector3f& v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline Vector3f operator+(const Vector3f& u, const Vector3f& v) {
    return Vector3f(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline Vector3f operator-(const Vector3f& u, const Vector3f& v) {
    return Vector3f(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline Vector3f operator*(const Vector3f& u, const Vector3f& v) {
    return Vector3f(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline Vector3f operator*(double t, const Vector3f& v) {
    return Vector3f(t * v.e[0], t * v.e[1], t * v.e[2]);
}

inline Vector3f operator*(const Vector3f& v, double t) {
    return t * v;
}

inline Vector3f operator/(Vector3f v, double t) {
    return (1 / t) * v;
}

inline double dot(const Vector3f& u, const Vector3f& v) {
    return u.e[0] * v.e[0]
        + u.e[1] * v.e[1]
        + u.e[2] * v.e[2];
}

inline Vector3f cross(const Vector3f& u, const Vector3f& v) {
    return Vector3f(u.e[1] * v.e[2] - u.e[2] * v.e[1],
        u.e[2] * v.e[0] - u.e[0] * v.e[2],
        u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline Vector3f unit_vector(Vector3f v) {
    return v / v.length();
}

inline Vector3f random_in_unit_sphere() {
    while (true) {
        auto p = Vector3f::random(-1, 1);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}

inline Vector3f random_unit_vector() {
    return unit_vector(random_in_unit_sphere());
}

inline Vector3f random_in_hemisphere(const Vector3f& normal) {
    Vector3f in_unit_sphere = random_in_unit_sphere();
    if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

inline Vector3f reflect(const Vector3f& v, const Vector3f& n) {
    return v - 2 * dot(v, n) * n;
}

inline Vector3f refract(const Vector3f& uv, const Vector3f& n, double etai_over_etat) {
    auto cos_theta = fmin(dot(-uv, n), 1.0);
    Vector3f r_out_perp = etai_over_etat * (uv + cos_theta * n);
    Vector3f r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
    return r_out_perp + r_out_parallel;
}

inline Vector3f random_in_unit_disk() {
    while (true) {
        auto p = Vector3f(random_double(-1, 1), random_double(-1, 1), 0);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}

inline Vector3f random_cosine_direction() {
    auto r1 = random_double();
    auto r2 = random_double();
    auto z = sqrt(1 - r2);

    auto phi = 2 * pi * r1;
    auto x = cos(phi) * sqrt(r2);
    auto y = sin(phi) * sqrt(r2);

    return Vector3f(x, y, z);
}