
#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>
#include "rtweekend.h"

using std::sqrt;

class Normal;
class point3;

class vec3 {
public:
    vec3() { x = y = z = 0; }
    vec3(Float e0, Float e1, Float e2) : x(e0), y(e1), z(e2) {}
    explicit vec3(const Normal& n);
    explicit vec3(const point3& p);
    vec3 Abs() { x = abs(x); y = abs(y); z = abs(z); }

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

    Float Length() const {
        return sqrt(LengthSquared());
    }

    Float LengthSquared() const {
        return x * x + y * y + z * z;
    }

    vec3 Normalize() {
        Float norm = Length();
        x /= norm;
        y /= norm;
        z /= norm;
        return (*this);
    }

    inline static vec3 random() {
        return vec3(RandomFloat(), RandomFloat(), RandomFloat());
    }

    inline static vec3 random(Float min, Float max) {
        return vec3(RandomFloat(min, max), RandomFloat(min, max), RandomFloat(min, max));
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
//using point3 = vec3;   // 3D point
using color = vec3;    // RGB color
//using Normal = vec3;
// vec3 Utility Functions

inline std::ostream& operator<<(std::ostream& out, const vec3& v) {
    return out << v.x << ' ' << v.y << ' ' << v.z;
}


inline vec3 Abs(vec3 u) {
    return vec3(abs(u.x), abs(u.y), abs(u.z));
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

inline vec3 Normalize(vec3 v) {
    return v / v.Length();
}

inline Float Dot(const vec3& u, const vec3& v) {
    return u.x * v.x
        + u.y * v.y
        + u.z * v.z;
}

inline vec3 Cross(const vec3& u, const vec3& v) {
    return vec3(u.y * v.z - u.z * v.y,
        u.z * v.x - u.x * v.z,
        u.x * v.y - u.y * v.x);
}

inline Float Distance(const vec3& u, const vec3& v)
{
    return (u - v).Length();
}

/*按xyz交换v的位置*/
inline vec3 Permute(const vec3& v, int x, int y, int z) {
    return vec3(v[x], v[y], v[z]);
}

inline Float MaxComponent(const vec3& v) {
    return std::max(std::max(v.x, v.y), v.z);
}

inline int MaxDimension(const vec3& v) {
    return (v.x > v.y) ? ((v.x > v.z) ? 0 : 2) : ((v.y > v.z) ? 1 : 2);
}

inline vec3 unit_vector(vec3 v) {
    return v / v.Length();
}


inline vec3 random_in_unit_sphere() {
        while (true) {
            auto p = vec3::random(-1, 1);
            if (p.LengthSquared() >= 1) continue;
            return p;
        }
    }

inline vec3 random_unit_vector() {
     return unit_vector(random_in_unit_sphere());
 }

 inline vec3 random_in_hemisphere(const vec3& normal) {
     vec3 in_unit_sphere = random_in_unit_sphere();
     if (Dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
         return in_unit_sphere;
     else
         return -in_unit_sphere;
 }
 
 inline vec3 reflect(const vec3& v, const vec3& n) {
     return v - 2 * Dot(v, n) * n;
 }

 inline vec3 refract(const vec3& uv, const vec3& n, Float etai_over_etat) {
     auto cos_theta = fmin(Dot(-uv, n), 1.0);
     vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
     vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.LengthSquared())) * n;
     return r_out_perp + r_out_parallel;
 }

 inline vec3 random_in_unit_disk() {
     while (true) {
         auto p = vec3(RandomFloat(-1, 1), RandomFloat(-1, 1), 0);
         if (p.LengthSquared() >= 1) continue;
         return p;
     }
 }

 inline vec3 RandomCosineDirection() {
     auto r1 = RandomFloat();
     auto r2 = RandomFloat();
     auto z = sqrt(1 - r2);

     auto phi = 2 * Pi * r1;
     auto x = cos(phi) * sqrt(r2);
     auto y = sin(phi) * sqrt(r2);

     return vec3(x, y, z);
 }

 inline void CoordinateSystem(const vec3& v1, vec3* v2, vec3* v3) {
     if (std::abs(v1.x) > std::abs(v1.y))
         *v2 = vec3(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
     else
         *v2 = vec3(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);
     *v3 = Cross(v1, *v2);
 }

 class point2 {
 public:
     point2() { x = y = 0; }
     point2(Float e0,Float e1):x(e0),y(e1){}
     point2 operator+(const point2& p) {
         return point2(x + p.x, y + p.y);
     }

     point2 operator*(Float a) {
         return point2(a * x, a * y);
     }
     Float x; Float y;
 };

 inline point2 operator*(Float a, point2 p) {
     return point2(a * p.x, a * p.y);
 }

 class Point2i {
 public:
     Point2i() { x = y = 0; };
     Point2i(int e1,int e2):x(e1),y(e2){}
     int x; int y;
 };

 class vec2 {
 public:
     vec2() { x = y = 0; }
     vec2(Float e0,Float e1):x(e0),y(e1){}

     Float& operator[](int i) {
         IN_RANGE((i == 0 || i == 1));
         if (i == 0) return this->x;
         if (i == 1) return this->y;
     }

     Float operator[](int i) const {
         IN_RANGE((i == 0 || i == 1));
         if (i == 0) return this->x;
         if (i == 1) return this->y;
     }

     Float x; Float y;
 };

 inline vec2 operator-(point2 p1, point2 p2) {
     return vec2(p1.x - p2.x, p1.y - p2.y);
 }

 class point3 {
 public:
     point3() { x = y = z = 0; }
     point3(Float e0, Float e1, Float e2) : x(e0), y(e1), z(e2) {}
     point3(vec3 v):x(v.x),y(v.y),z(v.z){}
     point3 Abs() { x = ::Abs(x); y = ::Abs(y); z = ::Abs(z); }

     point3 operator-() const { return point3(-x, -y, -z); }
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

     point3 operator+(const vec3& v) {
         return point3(x + v.x, y + v.y, z + v.z);
     }

     point3& operator+=(const vec3& v) {
         x += v.x;
         y += v.y;
         z += v.z;
         return *this;
     }
     point3 operator-(const vec3& v) {
         return point3(x - v.x, y - v.y, z - v.z);
     }

     point3& operator-=(const vec3& v)
     {
         return *this += (-v);
     }

     point3& operator*=(Float f)
     {
         x *= f; y *= f; z *= f;
         return *this;
     }

     Float Length() const {
         return sqrt(LengthSquared());
     }

     Float LengthSquared() const {
         return x * x + y * y + z * z;
     }

     inline static point3 random() {
         return point3(RandomFloat(), RandomFloat(), RandomFloat());
     }

     inline static point3 random(Float min, Float max) {
         return point3(RandomFloat(min, max), RandomFloat(min, max), RandomFloat(min, max));
     }

     bool near_zero() const {
         // Return true if the vector is close to zero in all dimensions.
         const auto s = 1e-8;
         return (fabs(x) < s) && (fabs(y) < s) && (fabs(z) < s);
     }

 public:
     Float x, y, z;
 };


 inline point3 Abs(point3 u) {
     return point3(Abs(u.x), Abs(u.y), Abs(u.z));
 }

 inline point3 operator+(const point3& u, const point3& v) {
     return point3(u.x + v.x, u.y + v.y, u.z + v.z);
 }

 inline vec3 operator-(const point3& u, const point3& v) {
     return vec3(u.x - v.x, u.y - v.y, u.z - v.z);
 }

 inline point3 operator+(const point3& u, const vec3& v) {
     return point3(u.x + v.x, u.y + v.y, u.z + v.z);
 }

 inline point3 operator-(const point3& u, const vec3& v) {
     return point3(u.x - v.x, u.y - v.y, u.z - v.z);
 }


 inline point3 operator*(const point3& u, const point3& v) {
     return point3(u.x * v.x, u.y * v.y, u.z * v.z);
 }

 inline point3 operator*(Float t, const point3& v) {
     return point3(t * v.x, t * v.y, t * v.z);
 }

 inline point3 operator*(const point3& v, Float t) {
     return t * v;
 }

 inline point3 operator/(point3 v, Float t) {
     return (1 / t) * v;
 }

 inline Float Distance(const point3& u, const point3& v)
 {
     return (u - v).Length();
 }

 inline Float DistanceSquared(const point3& u, const point3& v)
 {
     return (u - v).LengthSquared();
 }

 class Normal {
 public:
     Normal() { x = y = z = 0; }
     Normal(Float e0, Float e1, Float e2) : x(e0), y(e1), z(e2) {}
     explicit Normal(const vec3& v) :x(v.x), y(v.y), z(v.z) {}
     Normal Abs() { x = ::Abs(x); y = ::Abs(y); z = ::Abs(z); }
     bool operator==(const Normal& n) const {
         return x == n.x && y == n.y && z == n.z;
     }
     bool operator!=(const Normal& n) const {
         return x != n.x || y != n.y || z != n.z;
     }

     Normal operator-() const { return Normal(-x, -y, -z); }
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

     Normal operator+(const Normal& n) {
         return Normal(x + n.x, y + n.y, z + n.z);
     }

     Normal& operator+=(const Normal& v) {
         x += v.x;
         y += v.y;
         z += v.z;
         return *this;
     }

     Normal operator-(const Normal& n) {
         return Normal(x - n.x, y - n.y, z - n.z);
     }

     Normal& operator-=(const Normal& v)
     {
         return *this += (-v);
     }

     Normal& operator*=(const Float t) {
         x *= t;
         y *= t;
         z *= t;
         return *this;
     }

     Normal& operator/=(const Float t) {
         ZERO_DENOMINATOR(t);
         return *this *= 1 / t;
     }
     Normal operator*(Float t) {
         return Normal(x * t, y * t, z * t);
     }
     Normal operator/(Float t) {
         return Normal(x / t, y / t, z / t);
     }

     Float Length() const {
         return sqrt(LengthSquared());
     }

     Float LengthSquared() const {
         return x * x + y * y + z * z;
     }

     inline static Normal random() {
         return Normal(RandomFloat(), RandomFloat(), RandomFloat());
     }

     inline static Normal random(Float min, Float max) {
         return Normal(RandomFloat(min, max), RandomFloat(min, max), RandomFloat(min, max));
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

 inline Normal operator*(Float a, Normal n) {
     return Normal(a * n.x, a * n.y, a * n.z);
 }

 inline vec3::vec3(const Normal& n) :x(n.x), y(n.y), z(n.z) {};
 inline vec3::vec3(const point3& p) :x(p.x), y(p.y), z(p.z) {};


 inline Normal Normalize(Normal n) { return n / n.Length(); }

 inline Float AbsDot(Normal n1, vec3 v2) {
     return std::abs(n1.x * v2.x + n1.y * v2.y + n1.z * v2.z);
 }

 inline Normal Faceforward(const Normal& n, const vec3& v)
 {
     return Dot(vec3(n), v) >= 0.f ? n : -n;
 }

 using Vector3f = vec3;
 using Point3f = point3;
 using Point2f = point2;
 using Normal3f = Normal;
 using Vector2f = vec2;

#endif