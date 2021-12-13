#pragma once
#include <cmath>
#include <iostream>
#include "random.h"
#include "check.h"
#include <assert.h>
const double PI = 3.1415926535897932385;
using std::sqrt;


/*class vec3 {
public:
    vec3() { x = 0; y = 0; z = 0; }
    vec3(double e0, double e1, double e2) { x = e0; y = e1; z = e2; }

    

    vec3 operator-() const { return vec3(-x, -y, -z); }
    double operator[](int i) const { 
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    double& operator[](int i) {
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }

    vec3& operator+=(const vec3& v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    vec3& operator*=(const double t) {
        x *= t;
        y *= t;
        z *= t;
        return *this;
    }

    vec3& operator/=(const double t) {
        return *this *= 1 / t;
    }

    double Length() const {
        return sqrt(LengthSquared());
    }

    double LengthSquared() const {
        return x * x + y * y + z * z;
    }

    inline static vec3 random() {
        return vec3(rand() / (RAND_MAX + 1.0), rand() / (RAND_MAX + 1.0), rand() / (RAND_MAX + 1.0));
    }

    inline static vec3 random(double min, double max) {
        return vec3(min + (max - min) * rand() / (RAND_MAX + 1.0),
            min + (max - min) * rand() / (RAND_MAX + 1.0),
            min + (max - min) * rand() / (RAND_MAX + 1.0));
    }

    bool near_zero() const {
        // Return true if the vector is close to zero in all dimensions.
        const auto s = 1e-8;
        return (fabs(x) < s) && (fabs(y) < s) && (fabs(z) < s);
    }

    double* ToArray() {
        double rgb[3] = { x,y,z };
        return rgb;
    }


public:
    double x, y, z;
};

// Type aliases for vec3
   // RGB color

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

inline vec3 operator*(double t, const vec3& v) {
    return vec3(t * v.x, t * v.y, t * v.z);
}

inline vec3 operator*(const vec3& v, double t) {
    return t * v;
}

inline vec3 operator/(vec3 v, double t) {
    return (1 / t) * v;
}

inline double Dot(const vec3& u, const vec3& v) {
    return u.x * v.x
        + u.y * v.y
        + u.z * v.z;
}

inline vec3 Cross(const vec3& u, const vec3& v) {
    return vec3(u.y * v.z - u.z * v.y,
        u.z * v.x - u.x * v.z,
        u.x * v.y - u.y * v.x);
}

inline vec3 Normalize(vec3 v) {
    return v / v.Length();
}



*/






template <typename T>
class Vector2;
template <typename T>
class Vector3;
template <typename T>
class Point3;
template <typename T>
class Point2;


template <typename T>
inline bool isNaN(const T x) {
    return std::isnan(x);
}
template <>
inline bool isNaN(const int x) {
    return false;
}
/*......Vector......*/
template <typename T>
class Vector2
{
public:
    // Vector2 Public Methods
    Vector2() { x = y = 0; }
    Vector2(T xx, T yy) : x(xx), y(yy) { assert(!HasNaNs()); }
    Vector2(const Point2<T>& p)
    {
        x = p.x;
        y = p.y;
        assert(!HasNaNs());
    }
    bool HasNaNs() const { return isNaN(x) || isNaN(y); }

    Vector2(const Point3<T>& p)
        : x(p.x), y(p.y) {
        assert(!HasNaNs());
    }

    Vector2(const Vector2<T>& v) {
        assert(!v.HasNaNs());
        x = v.x;
        y = v.y;

    }
    Vector2<T>& operator=(const Vector2<T>& v) {
        assert(!v.HasNaNs());
        x = v.x;
        y = v.y;
        return *this;
    }

    Vector2<T> operator+(const Vector2<T>& v) const {
        assert(!v.HasNaNs());
        return Vector2(x + v.x, y + v.y);
    }

    Vector2<T>& operator+=(const Vector2<T>& v) {
        assert(!v.HasNaNs());
        x += v.x;
        y += v.y;
        return *this;
    }
    Vector2<T> operator-(const Vector2<T>& v) const {
        assert(!v.HasNaNs());
        return Vector2(x - v.x, y - v.y);
    }

    Vector2<T>& operator-=(const Vector2<T>& v) {
        assert(!v.HasNaNs());
        x -= v.x;
        y -= v.y;
        return *this;
    }
    bool operator==(const Vector2<T>& v) const { return x == v.x && y == v.y; }
    bool operator!=(const Vector2<T>& v) const { return x != v.x || y != v.y; }
    template <typename U>
    Vector2<T> operator*(U f) const {
        return Vector2<T>(f * x, f * y);
    }

    template <typename U>
    Vector2<T>& operator*=(U f) {
        assert(!isNaN(f));
        x *= f;
        y *= f;
        return *this;
    }
    template <typename U>
    Vector2<T> operator/(U f) const {
        assert(f != 0);
        double inv = (double)1 / f;
        return Vector2<T>(x * inv, y * inv);
    }

    template <typename U>
    Vector2<T>& operator/=(U f) {
        assert(f != 0);
        double inv = (double)1 / f;
        x *= inv;
        y *= inv;
        return *this;
    }
    Vector2<T> operator-() const { return Vector2<T>(-x, -y); }
    T operator[](int i) const {
        assert(i >= 0 && i <= 1);
        if (i == 0) return x;
        return y;
    }

    T& operator[](int i) {
        assert(i >= 0 && i <= 1);
        if (i == 0) return x;
        return y;
    }
    double LengthSquared() const { return x * x + y * y; }
    double Length() const { return sqrt(LengthSquared()); }


    T x, y;
};

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Vector2<T>& v) {
    os << "[ " << v.x << ", " << v.y << " ]";
    return os;
}

template <typename T>
class Vector3 {
public:
    // Vector3 Public Methods
    T operator[](int i) const {
        assert(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    T& operator[](int i) {
        assert(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    Vector3() { x = y = z = 0; }
    Vector3(T x, T y, T z) : x(x), y(y), z(z) { assert(!HasNaNs()); }
    bool HasNaNs() const { return isNaN(x) || isNaN(y) || isNaN(z); }
    /*Vector3(const Point3<T>& p) {
        assert(!v.HasNaNs());
        x = p.x;
        y = p.y;
        z = p.z;
    }*/
    explicit Vector3(const Point3<T>& p)
        : x(p.x), y(p.y), z(p.z) {
        assert(!HasNaNs());
    }

    template <typename U>
    inline Vector3<T> MutiXYZ(const Vector3<U>& v) {
        Vector3<T> tmp(x * v.x, y * v.y, z * v.z);
        return tmp;
    }

    Vector3(const Vector3<T>& v) {
        assert(!v.HasNaNs());
        x = v.x;
        y = v.y;
        z = v.z;
    }

    Vector3<T>& operator=(const Vector3<T>& v) {
        assert(!v.HasNaNs());
        x = v.x;
        y = v.y;
        z = v.z;
        return *this;
    }

    Vector3<T> operator+(const Vector3<T>& v) const {
        assert(!v.HasNaNs());
        return Vector3(x + v.x, y + v.y, z + v.z);
    }
    Vector3<T>& operator+=(const Vector3<T>& v) {
        assert(!v.HasNaNs());
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    Vector3<T> operator-(const Vector3<T>& v) const {
        assert(!v.HasNaNs());
        return Vector3(x - v.x, y - v.y, z - v.z);
    }
    Vector3<T>& operator-=(const Vector3<T>& v) {
        assert(!v.HasNaNs());
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
    bool operator==(const Vector3<T>& v) const {
        return x == v.x && y == v.y && z == v.z;
    }
    bool operator!=(const Vector3<T>& v) const {
        return x != v.x || y != v.y || z != v.z;
    }
    template <typename U>
    Vector3<T> operator*(U s) const {
        return Vector3<T>(s * x, s * y, s * z);
    }
    template <typename U>
    Vector3<T> XYZMuti(const Vector3<U>& s) const {
        return Vector3<T>(s.x * x, s.y * y, s.z * z);
    }

    template <typename U>
    Vector3<T>& operator*=(U s) {
        assert(!isNaN(s));
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }
    template <typename U>
    Vector3<T> operator/(U f) const {
        assert(f != 0);
        double inv = (double)1 / f;
        return Vector3<T>(x * inv, y * inv, z * inv);
    }

    template <typename U>
    Vector3<T>& operator/=(U f) {
        assert(f != 0);
        double inv = (double)1 / f;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }
    Vector3<T> operator-() const { return Vector3<T>(-x, -y, -z); }
    double LengthSquared() const { return x * x + y * y + z * z; }
    double Length() const { return sqrt(LengthSquared()); }


    inline static Vector3<T> random() {
        return Vector3<T>(rand() / (RAND_MAX + 1.0), rand() / (RAND_MAX + 1.0), rand() / (RAND_MAX + 1.0));
    }

    inline static Vector3<T> random(double min, double max) {
        return Vector3<T>(min + (max - min) * rand() / (RAND_MAX + 1.0),
            min + (max - min) * rand() / (RAND_MAX + 1.0),
            min + (max - min) * rand() / (RAND_MAX + 1.0));
    }

    bool near_zero() const {
        // Return true if the vector is close to zero in all dimensions.
        const auto s = 1e-8;
        return (fabs(x) < s) && (fabs(y) < s) && (fabs(z) < s);
    }

    // Vector3 Public Data

    T* ToArray() {
        T xyz[3] = { x,y,z };
        return xyz;
    }
    T x, y, z;
};
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Vector3<T>& v) {
    os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
    return os;
}

template <typename T, typename U>
inline Vector3<T> operator*(U t, const Vector3<T>& v) {
    return Vector3<T>(t * v.x, t * v.y, t * v.z);
}


template <typename T>
inline T Dot(const Vector3<T>& v1, const Vector3<T>& v2) {
    assert(!v1.HasNaNs() && !v2.HasNaNs());
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <typename T>
Vector3<T> Abs(const Vector3<T>& v) {
    return Vector3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

template <typename T>
inline T AbsDot(const Vector3<T>& v1, const Vector3<T>& v2) {
    assert(!v1.HasNaNs() && !v2.HasNaNs());
    return std::abs(Dot(v1, v2));
}

template <typename T>
inline Vector3<T> Cross(const Vector3<T>& v1, const Vector3<T>& v2) {
    assert(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
        (v1x * v2y) - (v1y * v2x));
}

template <typename T>
inline Vector3<T> Normalize(const Vector3<T>& v) {
    return v / v.Length();
}
template <typename T>
T MinComponent(const Vector3<T>& v) {
    return std::min(v.x, std::min(v.y, v.z));
}

template <typename T>
T MaxComponent(const Vector3<T>& v) {
    return std::max(v.x, std::max(v.y, v.z));
}

template <typename T>
int MaxDimension(const Vector3<T>& v) {
    return (v.x > v.y) ? ((v.x > v.z) ? 0 : 2) : ((v.y > v.z) ? 1 : 2);
}

template <typename T>
Vector3<T> Min(const Vector3<T>& p1, const Vector3<T>& p2) {
    return Vector3<T>(std::min(p1.x, p2.x), std::min(p1.y, p2.y),
        std::min(p1.z, p2.z));
}

template <typename T>
Vector3<T> Max(const Vector3<T>& p1, const Vector3<T>& p2) {
    return Vector3<T>(std::max(p1.x, p2.x), std::max(p1.y, p2.y),
        std::max(p1.z, p2.z));
}

template <typename T>
Vector3<T> Permute(const Vector3<T>& v, int x, int y, int z) {
    return Vector3<T>(v[x], v[y], v[z]);
}

template <typename T>
inline void CoordinateSystem(const Vector3<T>& v1, Vector3<T>* v2,
    Vector3<T>* v3) {
    if (std::abs(v1.x) > abs(v1.y))
        *v2 = Vector3<T>(-v1.z, 0, v1.x) / sqrt(v1.x * v1.x + v1.z * v1.z);
    else
        *v2 = Vector3<T>(0, v1.z, -v1.y) / sqrt(v1.y * v1.y + v1.z * v1.z);
    *v3 = Cross(v1, *v2);
}

template <typename T, typename U>
inline Vector2<T> operator*(U f, const Vector2<T>& v) {
    return v * f;
}
template <typename T>
inline double Dot(const Vector2<T>& v1, const Vector2<T>& v2) {
    assert(!v1.HasNaNs() && !v2.HasNaNs());
    return v1.x * v2.x + v1.y * v2.y;
}

template <typename T>
inline double AbsDot(const Vector2<T>& v1, const Vector2<T>& v2) {
    assert(!v1.HasNaNs() && !v2.HasNaNs());
    return std::abs(Dot(v1, v2));
}

template <typename T>
inline Vector2<T> Normalize(const Vector2<T>& v) {
    return v / v.Length();
}
template <typename T>
Vector2<T> Abs(const Vector2<T>& v) {
    return Vector2<T>(std::abs(v.x), std::abs(v.y));
}


typedef Vector2<double> Vector2f;
typedef Vector2<int> Vector2i;
typedef Vector3<double> Vector3f;
typedef Vector3<int> Vector3i;


template <typename T>
class Point2 {
public:
    // Point2 Public Methods
    explicit Point2(const Point3<T>& p) {
        x = p.x;
        y = p.y;
        assert(!HasNaNs());
    }
    Point2() { x = y = 0; }
    Point2(T xx, T yy) : x(xx), y(yy) { assert(!HasNaNs()); }

    template <typename U>
    explicit Point2(const Point2<U>& p) {
        x = (T)p.x;
        y = (T)p.y;
        assert(!HasNaNs());
    }

    template <typename U>
    explicit Point2(const Vector2<U>& p) {
        x = (T)p.x;
        y = (T)p.y;
        assert(!HasNaNs());
    }

    template <typename U>
    explicit operator Vector2<U>() const {
        return Vector2<U>(x, y);
    }


    Point2(const Point2<T>& p) {
        assert(!p.HasNaNs());
        x = p.x;
        y = p.y;
    }

    Point2<T>& operator=(const Point2<T>& p) {
        assert(!p.HasNaNs());
        x = p.x;
        y = p.y;
        return *this;
    }

    Point2<T> operator+(const Vector2<T>& v) const {
        assert(!v.HasNaNs());
        return Point2<T>(x + v.x, y + v.y);
    }

    Point2<T>& operator+=(const Vector2<T>& v) {
        assert(!v.HasNaNs());
        x += v.x;
        y += v.y;
        return *this;
    }
    Vector2<T> operator-(const Point2<T>& p) const {
        assert(!p.HasNaNs());
        return Vector2<T>(x - p.x, y - p.y);
    }

    Point2<T> operator-(const Vector2<T>& v) const {
        assert(!v.HasNaNs());
        return Point2<T>(x - v.x, y - v.y);
    }
    Point2<T> operator-() const { return Point2<T>(-x, -y); }
    Point2<T>& operator-=(const Vector2<T>& v) {
        assert(!v.HasNaNs());
        x -= v.x;
        y -= v.y;
        return *this;
    }
    Point2<T>& operator+=(const Point2<T>& p) {
        assert(!p.HasNaNs());
        x += p.x;
        y += p.y;
        return *this;
    }
    Point2<T> operator+(const Point2<T>& p) const {
        assert(!p.HasNaNs());
        return Point2<T>(x + p.x, y + p.y);
    }
    template <typename U>
    Point2<T> operator*(U f) const {
        return Point2<T>(f * x, f * y);
    }
    template <typename U>
    Point2<T>& operator*=(U f) {
        x *= f;
        y *= f;
        return *this;
    }
    template <typename U>
    Point2<T> operator/(U f) const {
        assert(f != 0);
        double inv = (double)1 / f;
        return Point2<T>(inv * x, inv * y);
    }
    template <typename U>
    Point2<T>& operator/=(U f) {
        assert(f != 0);
        double inv = (double)1 / f;
        x *= inv;
        y *= inv;
        return *this;
    }
    T operator[](int i) const {
        assert(i >= 0 && i <= 1);
        if (i == 0) return x;
        return y;
    }

    T& operator[](int i) {
        assert(i >= 0 && i <= 1);
        if (i == 0) return x;
        return y;
    }
    bool operator==(const Point2<T>& p) const { return x == p.x && y == p.y; }
    bool operator!=(const Point2<T>& p) const { return x != p.x || y != p.y; }
    bool HasNaNs() const { return isNaN(x) || isNaN(y); }

    // Point2 Public Data
    T x, y;
};

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Point2<T>& v) {
    os << "[ " << v.x << ", " << v.y << " ]";
    return os;
}

template <typename T>
class Point3 {
public:
    // Point3f Public Methods
    Point3() { x = y = z = 0; }
    Point3(T x, T y, T z) : x(x), y(y), z(z) { assert(!HasNaNs()); }
    bool HasNaNs() const { return isNaN(x) || isNaN(y) || isNaN(z); }
    template <typename U>
    Point3(const Point3<U>& p)
        : x((T)p.x), y((T)p.y), z((T)p.z) {
        assert(!HasNaNs());
    }
    
    template <typename U>
    operator Vector3<U>() const {
        return Vector3<U>(x, y, z);
    }

    explicit Point3(const Vector3<T>& p) {
        assert(!p.HasNaNs());
        x = p.x;
        y = p.y;
        z = p.z;
    }

    Point3(const Point3<T>& p) {
        assert(!p.HasNaNs());
        x = p.x;
        y = p.y;
        z = p.z;
    }

    Point3<T>& operator=(const Point3<T>& p) {
        assert(!p.HasNaNs());
        x = p.x;
        y = p.y;
        z = p.z;
        return *this;
    }

    Point3<T> operator+(const Vector3<T>& v) const {
        assert(!v.HasNaNs());
        return Point3<T>(x + v.x, y + v.y, z + v.z);
    }
    Point3<T>& operator+=(const Vector3<T>& v) {
        assert(!v.HasNaNs());
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    /*Vector3<T> operator-(const Point3<T>& p) const {
        assert(!p.HasNaNs());
        return Vector3<T>(x - p.x, y - p.y, z - p.z);
    }*/
    Point3<T> operator-(const Vector3<T>& v) const {
        assert(!v.HasNaNs());
        return Point3<T>(x - v.x, y - v.y, z - v.z);
    }
    Point3<T>& operator-=(const Vector3<T>& v) {
        assert(!v.HasNaNs());
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
    Point3<T>& operator+=(const Point3<T>& p) {
        assert(!p.HasNaNs());
        x += p.x;
        y += p.y;
        z += p.z;
        return *this;
    }
    Point3<T> operator+(const Point3<T>& p) const {
        assert(!p.HasNaNs());
        return Point3<T>(x + p.x, y + p.y, z + p.z);
    }
    Point3<T> operator-(const Point3<T>& p) const {
        assert(!p.HasNaNs());
        return Point3<T>(x - p.x, y - p.y, z - p.z);
    } 

    /*Point3<T>& operator-=(const Point3<T>& v) {
        assert(!v.HasNaNs());
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }*/
    template <typename U>
    Point3<T> operator*(U f) const {
        return Point3<T>(f * x, f * y, f * z);
    }
    template <typename U>
    Point3<T>& operator*=(U f) {
        x *= f;
        y *= f;
        z *= f;
        return *this;
    }
    template <typename U>
    Point3<T> operator/(U f) const {
        assert(f != 0);
        double inv = (double)1 / f;
        return Point3<T>(inv * x, inv * y, inv * z);
    }
    template <typename U>
    Point3<T>& operator/=(U f) {
        assert(f != 0);
        double inv = (double)1 / f;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }
    T operator[](int i) const {
        assert(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }

    T& operator[](int i) {
        assert(i >= 0 && i <= 2);
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    bool operator==(const Point3<T>& p) const {
        return x == p.x && y == p.y && z == p.z;
    }
    bool operator!=(const Point3<T>& p) const {
        return x != p.x || y != p.y || z != p.z;
    }

    Point3<T> operator-() const { return Point3<T>(-x, -y, -z); }

    double LengthSquared() const { return x * x + y * y + z * z; }
    double Length() const { return sqrt(LengthSquared()); }


    inline static Point3<T> random() {
        return Point3<T>(rand() / (RAND_MAX + 1.0), rand() / (RAND_MAX + 1.0), rand() / (RAND_MAX + 1.0));
    }

    inline static Point3<T> random(double min, double max) {
        return Point3<T>(min + (max - min) * rand() / (RAND_MAX + 1.0),
            min + (max - min) * rand() / (RAND_MAX + 1.0),
            min + (max - min) * rand() / (RAND_MAX + 1.0));
    }

    bool near_zero() const {
        // Return true if the vector is close to zero in all dimensions.
        const auto s = 1e-8;
        return (fabs(x) < s) && (fabs(y) < s) && (fabs(z) < s);
    }

    // Point3f Public Data
    
    T* ToArray() {
        T xyz[3] = { x,y,z };
        return xyz;
    }
    
    T x, y, z;
};

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Point3<T>& v) {
    os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
    return os;
}

typedef Point2<double> Point2f;
typedef Point2<int> Point2i;
typedef Point3<double> Point3d;
typedef Point3<int> Point3i;


template <typename T, typename U>
inline Point3<T> operator*(U f, const Point3<T>& p) {
    assert(!p.HasNaNs());
    return p * f;
}

template <typename T>
Point3<T> Lerp(double t, const Point3<T>& p0, const Point3<T>& p1) {
    return (1 - t) * p0 + t * p1;
}

template <typename T>
Point3<T> Min(const Point3<T>& p1, const Point3<T>& p2) {
    return Point3<T>(std::min(p1.x, p2.x), std::min(p1.y, p2.y),
        std::min(p1.z, p2.z));
}

template <typename T>
Point3<T> Max(const Point3<T>& p1, const Point3<T>& p2) {
    return Point3<T>(std::max(p1.x, p2.x), std::max(p1.y, p2.y),
        std::max(p1.z, p2.z));
}

template <typename T>
Point3<T> Floor(const Point3<T>& p) {
    return Point3<T>(std::floor(p.x), std::floor(p.y), std::floor(p.z));
}

template <typename T>
Point3<T> Ceil(const Point3<T>& p) {
    return Point3<T>(std::ceil(p.x), std::ceil(p.y), std::ceil(p.z));
}

template <typename T>
Point3<T> Abs(const Point3<T>& p) {
    return Point3<T>(std::abs(p.x), std::abs(p.y), std::abs(p.z));
}

template <typename T>
inline double Length(const Point2<T>& p1, const Point2<T>& p2) {
    return (p1 - p2).Length();
}

template <typename T>
inline double LengthSquared(const Point2<T>& p1, const Point2<T>& p2) {
    return (p1 - p2).LengthSquared();
}

template <typename T, typename U>
inline Point2<T> operator*(U f, const Point2<T>& p) {
    assert(!p.HasNaNs());
    return p * f;
}

template <typename T>
Point2<T> Floor(const Point2<T>& p) {
    return Point2<T>(std::floor(p.x), std::floor(p.y));
}

template <typename T>
Point2<T> Ceil(const Point2<T>& p) {
    return Point2<T>(std::ceil(p.x), std::ceil(p.y));
}

template <typename T>
Point2<T> Lerp(double t, const Point2<T>& v0, const Point2<T>& v1) {
    return (1 - t) * v0 + t * v1;
}

template <typename T>
Point2<T> Min(const Point2<T>& pa, const Point2<T>& pb) {
    return Point2<T>(std::min(pa.x, pb.x), std::min(pa.y, pb.y));
}

template <typename T>
Point2<T> Max(const Point2<T>& pa, const Point2<T>& pb) {
    return Point2<T>(std::max(pa.x, pb.x), std::max(pa.y, pb.y));
}

template <typename T>
Point3<T> Permute(const Point3<T>& p, int x, int y, int z) {
    return Point3<T>(p[x], p[y], p[z]);
}
template <typename T>
inline double Length(const Point3<T>& p1, const Point3<T>& p2) {
    return (p1 - p2).Length();
}

template <typename T>
inline double LengthSquared(const Point3<T>& p1, const Point3<T>& p2) {
    return (p1 - p2).LengthSquared();
}

template <typename T>
inline double Dot(const Point3<T>& v1, const Point3<T>& v2) {
    assert(!v1.HasNaNs() && !v2.HasNaNs());
    return v1.x * v2.x + v1.y * v2.y;
}

template <typename T>
inline double AbsDot(const Point3<T>& v1, const Point3<T>& v2) {
    assert(!v1.HasNaNs() && !v2.HasNaNs());
    return std::abs(Dot(v1, v2));
}

template <typename T>
inline Point3<T> Normalize(const Point3<T>& v) {
    return v / v.Length();
}
template <typename T>
inline Point3<T> Cross(const Point3<T>& v1, const Point3<T>& v2) {
    assert(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
        (v1x * v2y) - (v1y * v2x));
}


inline Vector3f random_in_unit_sphere() {
    while (true) {
        auto p = Vector3f::random(-1, 1);
        if (p.LengthSquared() >= 1) continue;
        return p;
    }
}

inline Vector3f random_unit_vector() {
    return Normalize(random_in_unit_sphere());
}

inline Vector3f random_in_hemisphere(const Vector3f& normal) {
    Vector3f in_unit_sphere = random_in_unit_sphere();
    if (Dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

inline Vector3f reflect(const Vector3f& v, const Vector3f& n) {
    return v - 2 * Dot(v, n) * n;
}

inline Vector3f refract(const Vector3f& uv, const Vector3f& n, double etai_over_etat) {
    auto cos_theta = fmin(Dot(-uv, n), 1.0);
    Vector3f r_out_perp = etai_over_etat * (uv + cos_theta * n);
    Vector3f r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.LengthSquared())) * n;
    return r_out_perp + r_out_parallel;
}

inline Vector3f random_in_unit_disk() {
    while (true) {
        auto p = Vector3f(random_double(-1, 1), random_double(-1, 1), 0);
        if (p.LengthSquared() >= 1) continue;
        return p;
    }
}

inline Vector3f random_cosine_direction() {
    auto r1 = random_double();
    auto r2 = random_double();
    auto z = sqrt(1 - r2);

    auto phi = 2 * PI * r1;
    auto x = cos(phi) * sqrt(r2);
    auto y = sin(phi) * sqrt(r2);

    return Vector3f(x, y, z);
}

//template <typename T>
inline int MaxDimension(const Vector3f& v) {
    return (v.x > v.y) ? ((v.x > v.z) ? 0 : 2) : ((v.y > v.z) ? 1 : 2);
}



using Point3f = Vector3f;   // 3D point
//using color = Vector3f;
using color = Vector3f;
using Normal3f = Vector3f;


