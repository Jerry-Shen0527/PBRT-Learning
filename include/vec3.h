#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>
#include<assert.h>
#include "rtweekend.h"
using std::sqrt;

template <typename T>
class Point2;
template <typename T>
class Point3;
template <typename T>
class Normal3;


template <typename T>
class Vector2 {
public:
    // Vector2 Public Methods
    Vector2() { x = y = 0; }
    Vector2(T xx, T yy) : x(xx), y(yy) { //DCHECK(!HasNaNs()); 
    }
    bool HasNaNs() const { return isNaN(x) || isNaN(y); }
    explicit Vector2(const Point2<T>& p) : x(p.x), y(p.y) {
        //DCHECK(!HasNaNs());
    }
    explicit Vector2(const Point3<T>& p) : x(p.x), y(p.y) {
        //DCHECK(!HasNaNs());
    }
#ifndef NDEBUG
    // The default versions of these are fine for release builds; for debug
    // we define them so that we can add the Assert checks.
    Vector2(const Vector2<T>& v) {
        DCHECK(!v.HasNaNs());
        x = v.x;
        y = v.y;
    }
    Vector2<T>& operator=(const Vector2<T>& v) {
        DCHECK(!v.HasNaNs());
        x = v.x;
        y = v.y;
        return *this;
    }
#endif  // !NDEBUG

    Vector2<T> operator+(const Vector2<T>& v) const {
        //DCHECK(!v.HasNaNs());
        return Vector2(x + v.x, y + v.y);
    }

    Vector2<T>& operator+=(const Vector2<T>& v) {
        //DCHECK(!v.HasNaNs());
        x += v.x;
        y += v.y;
        return *this;
    }
    Vector2<T> operator-(const Vector2<T>& v) const {
        //DCHECK(!v.HasNaNs());
        return Vector2(x - v.x, y - v.y);
    }

    Vector2<T>& operator-=(const Vector2<T>& v) {
        //DCHECK(!v.HasNaNs());
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
        //DCHECK(!isNaN(f));
        x *= f;
        y *= f;
        return *this;
    }
    template <typename U>
    Vector2<T> operator/(U f) const {
        //CHECK_NE(f, 0);
        assert(f != 0);
        Float inv = (Float)1 / f;
        return Vector2<T>(x * inv, y * inv);
    }

    template <typename U>
    Vector2<T>& operator/=(U f) {
        //CHECK_NE(f, 0);
        assert(f != 0);
        Float inv = (Float)1 / f;
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
    Float LengthSquared() const { return x * x + y * y; }
    Float Length() const { return std::sqrt(LengthSquared()); }

    // Vector2 Public Data
    T x, y;
};

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Vector2<T>& v) {
    os << "[ " << v.x << ", " << v.y << " ]";
    return os;
}

template <>
inline std::ostream& operator<<(std::ostream& os, const Vector2<Float>& v) {
    //os << StringPrintf("[ %f, %f ]", v.x, v.y);
    os << "[ " << v.x << ", " << v.y << " ]";
    return os;
}


template <typename T>
class Vector3 {
public:
    Vector3() { x = y = z = 0; }
    Vector3(T x, T y, T z) : x(x), y(y), z(z) { //DCHECK(!HasNaNs()); 
    }
    bool HasNaNs() const { return isNaN(x) || isNaN(y) || isNaN(z); }
    explicit Vector3(const Point3<T>& p) : x(p.x), y(p.y), z(p.z) {
        DCHECK(!HasNaNs());
    }
    explicit Vector3(const Normal3<T>& n) : x(n.x), y(n.y), z(n.z) {
        //DCHECK(!n.HasNaNs());
    }

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
    

#ifndef NDEBUG
    // The default versions of these are fine for release builds; for debug
    // we define them so that we can add the Assert checks.
    Vector3(const Vector3<T>& v) {
        DCHECK(!v.HasNaNs());
        x = v.x;
        y = v.y;
        z = v.z;
    }

    Vector3<T>& operator=(const Vector3<T>& v) {
        DCHECK(!v.HasNaNs());
        x = v.x;
        y = v.y;
        z = v.z;
        return *this;
    }
#endif  // !NDEBUG
    Vector3<T> operator+(const Vector3<T>& v) const {
        //DCHECK(!v.HasNaNs());
        return Vector3(x + v.x, y + v.y, z + v.z);
    }
    Vector3<T>& operator+=(const Vector3<T>& v) {
        //DCHECK(!v.HasNaNs());
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    Vector3<T> operator-(const Vector3<T>& v) const {
        //DCHECK(!v.HasNaNs());
        return Vector3(x - v.x, y - v.y, z - v.z);
    }
    Vector3<T>& operator-=(const Vector3<T>& v) {
        //DCHECK(!v.HasNaNs());
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
    //add auto
    template <typename U>
    Vector3<T> operator*(Vector3<U> s) const {
        return Vector3<T>(s.x * x, s.y * y, s.z * z);
    }

    template <typename U>
    Vector3<T> operator*(U s) const {
        return Vector3<T>(s * x, s * y, s * z);
    }
    template <typename U>
    Vector3<T>& operator*=(U s) {
        //DCHECK(!isNaN(s));
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }
    template <typename U>
    Vector3<T> operator/(U f) const {
        //CHECK_NE(f, 0);
        assert(f != 0);
        Float inv = (Float)1 / f;
        return Vector3<T>(x * inv, y * inv, z * inv);
    }

    template <typename U>
    Vector3<T>& operator/=(U f) {
        //CHECK_NE(f, 0);
        assert(f != 0);
        Float inv = (Float)1 / f;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }
    Vector3<T> operator-() const { return Vector3<T>(-x, -y, -z); }
    Float LengthSquared() const { return x * x + y * y + z * z; }
    Float Length() const { return std::sqrt(LengthSquared()); }
   

public:
    inline static Vector3<Float> random() {
        return Vector3(random_Float(), random_Float(), random_Float());
    }

    inline static Vector3<Float> random(double min, double max) {
        return Vector3(random_Float(min, max), random_Float(min, max), random_Float(min, max));
    }

    bool near_zero() const {
        // Return true if the vector is close to zero in all dimensions.
        const auto s = 1e-8;
        return (fabs(x) < s) && (fabs(y) < s) && (fabs(z) < s);
    }

    // Vector3 Public Data
    T x, y, z;
};

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Vector3<T>& v) {
    os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
    return os;
}


template <>
inline std::ostream& operator<<(std::ostream& os, const Vector3<Float>& v) {
    //os << StringPrintf("[ %f, %f, %f ]", v.x, v.y, v.z);
    os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
    return os;
}

// Type aliases for vec3
using Vector3f = Vector3<Float>;
using Vector2f = Vector2<Float>;
//using point3 = vec3;   // 3D point
using color = Vector3<Float>;    // RGB color


template <typename T>
class Point2 {
public:
    // Point2 Public Methods
    Point2() { x = y = 0; }
    Point2(T xx, T yy) : x(xx), y(yy) { //DCHECK(!HasNaNs()); 
    }

    explicit Point2(const Point3<T>& p) : x(p.x), y(p.y) { //DCHECK(!HasNaNs()); 
    }

    template <typename U>
    explicit Point2(const Point2<U>& p) {
        x = (T)p.x;
        y = (T)p.y;
        //DCHECK(!HasNaNs());
    }

    template <typename U>
    explicit Point2(const Vector2<U>& p) {
        x = (T)p.x;
        y = (T)p.y;
        //DCHECK(!HasNaNs());
    }

    template <typename U>
    explicit operator Vector2<U>() const {
        return Vector2<U>(x, y);
    }

#ifndef NDEBUG
    Point2(const Point2<T>& p) {
        DCHECK(!p.HasNaNs());
        x = p.x;
        y = p.y;
    }

    Point2<T>& operator=(const Point2<T>& p) {
        DCHECK(!p.HasNaNs());
        x = p.x;
        y = p.y;
        return *this;
    }
#endif  // !NDEBUG
    Point2<T> operator+(const Vector2<T>& v) const {
        //DCHECK(!v.HasNaNs());
        return Point2<T>(x + v.x, y + v.y);
    }

    Point2<T>& operator+=(const Vector2<T>& v) {
        //DCHECK(!v.HasNaNs());
        x += v.x;
        y += v.y;
        return *this;
    }
    Vector2<T> operator-(const Point2<T>& p) const {
        //DCHECK(!p.HasNaNs());
        return Vector2<T>(x - p.x, y - p.y);
    }

    Point2<T> operator-(const Vector2<T>& v) const {
        //DCHECK(!v.HasNaNs());
        return Point2<T>(x - v.x, y - v.y);
    }
    Point2<T> operator-() const { return Point2<T>(-x, -y); }
    Point2<T>& operator-=(const Vector2<T>& v) {
        //DCHECK(!v.HasNaNs());
        x -= v.x;
        y -= v.y;
        return *this;
    }
    Point2<T>& operator+=(const Point2<T>& p) {
        //DCHECK(!p.HasNaNs());
        x += p.x;
        y += p.y;
        return *this;
    }
    Point2<T> operator+(const Point2<T>& p) const {
        //DCHECK(!p.HasNaNs());
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
        //CHECK_NE(f, 0);
        assert(f != 0);
        Float inv = (Float)1 / f;
        return Point2<T>(inv * x, inv * y);
    }
    template <typename U>
    Point2<T>& operator/=(U f) {
        //CHECK_NE(f, 0);
        assert(f != 0);
        Float inv = (Float)1 / f;
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

template <>
inline std::ostream& operator<<(std::ostream& os, const Point2<Float>& v) {
    //os << StringPrintf("[ %f, %f ]", v.x, v.y);
    os << "[ " << v.x << ", " << v.y << " ]";
    return os;
}

template <typename T>
class Point3 {
public:
    // Point3 Public Methods
    Point3() { x = y = z = 0; }
    Point3(T x, T y, T z) : x(x), y(y), z(z) { //CHECK(!HasNaNs());
    }

    template <typename U>
    explicit Point3(const Point3<U>& p)
        : x((T)p.x), y((T)p.y), z((T)p.z) {
       // DCHECK(!HasNaNs());
    }
    template <typename U>
    explicit operator Vector3<U>() const {
        return Vector3<U>(x, y, z);
    }
#ifndef NDEBUG
    Point3(const Point3<T>& p) {
        DCHECK(!p.HasNaNs());
        x = p.x;
        y = p.y;
        z = p.z;
    }

    Point3<T>& operator=(const Point3<T>& p) {
        DCHECK(!p.HasNaNs());
        x = p.x;
        y = p.y;
        z = p.z;
        return *this;
    }
#endif  // !NDEBUG
    Point3<T> operator+(const Vector3<T>& v) const {
        //DCHECK(!v.HasNaNs());
        return Point3<T>(x + v.x, y + v.y, z + v.z);
    }
    Point3<T>& operator+=(const Vector3<T>& v) {
        //DCHECK(!v.HasNaNs());
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    Vector3<T> operator-(const Point3<T>& p) const {
        //DCHECK(!p.HasNaNs());
        return Vector3<T>(x - p.x, y - p.y, z - p.z);
    }
    Point3<T> operator-(const Vector3<T>& v) const {
        //DCHECK(!v.HasNaNs());
        return Point3<T>(x - v.x, y - v.y, z - v.z);
    }
    Point3<T>& operator-=(const Vector3<T>& v) {
        //DCHECK(!v.HasNaNs());
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
    Point3<T>& operator+=(const Point3<T>& p) {
        //DCHECK(!p.HasNaNs());
        x += p.x;
        y += p.y;
        z += p.z;
        return *this;
    }
    Point3<T> operator+(const Point3<T>& p) const {
        //DCHECK(!p.HasNaNs());
        return Point3<T>(x + p.x, y + p.y, z + p.z);
    }
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
        //CHECK_NE(f, 0);
        assert(f != 0);
        Float inv = (Float)1 / f;
        return Point3<T>(inv * x, inv * y, inv * z);
    }
    template <typename U>
    Point3<T>& operator/=(U f) {
        //CHECK_NE(f, 0);
        assert(f != 0);
        Float inv = (Float)1 / f;
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
    bool HasNaNs() const { return isNaN(x) || isNaN(y) || isNaN(z); }
    Point3<T> operator-() const { return Point3<T>(-x, -y, -z); }

    //add auto
    inline static Point3<Float> random() {
        return Point3(random_Float(), random_Float(), random_Float());
    }

    inline static Point3<Float> random(double min, double max) {
        return Point3(random_Float(min, max), random_Float(min, max), random_Float(min, max));
    }
    // Point3 Public Data
    T x, y, z;
};
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Point3<T>& v) {
    os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
    return os;
}


template <>
inline std::ostream& operator<<(std::ostream& os, const Point3<Float>& v) {
    //os << StringPrintf("[ %f, %f, %f ]", v.x, v.y, v.z);
    os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
    return os;
}

using Point3f = Point3<Float>;
using Point2f = Point2<Float>;

/*
class vec3 {
public:
    vec3() : e{ 0,0,0 } {}
    vec3(double e0, double e1, double e2) : e{ e0, e1, e2 } {}

    double x() const { return e[0]; }//const indicate this to be const which means invalid to modify.
    double y() const { return e[1]; }
    double z() const { return e[2]; }

    vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
    double operator[](int i) const { return e[i]; }
    double& operator[](int i) { return e[i]; }//could be left value.

    vec3& operator+=(const vec3& v) {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    vec3& operator*=(const double t) {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    vec3& operator/=(const double t) {
        return *this *= 1 / t;
    }

    double length() const {
        return sqrt(length_squared());
    }

    double length_squared() const {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }

public:
    //rand_vec3
    inline static vec3 random() {
        return vec3(random_double(), random_double(), random_double());
    }

    inline static vec3 random(double min, double max) {
        return vec3(random_double(min, max), random_double(min, max), random_double(min, max));
    }

public:
    bool near_zero() const {
        // Return true if the vector is close to zero in all dimensions.
        const auto s = 1e-8;
        return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
    }

public:
    double e[3];
};
*/
template <typename T>
class Normal3 {
public:
    // Normal3 Public Methods
    Normal3() { x = y = z = 0; }
    Normal3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) { //DCHECK(!HasNaNs()); 
    }
    Normal3<T> operator-() const { return Normal3(-x, -y, -z); }
    Normal3<T> operator+(const Normal3<T>& n) const {
        //DCHECK(!n.HasNaNs());
        return Normal3<T>(x + n.x, y + n.y, z + n.z);
    }

    Normal3<T>& operator+=(const Normal3<T>& n) {
        //DCHECK(!n.HasNaNs());
        x += n.x;
        y += n.y;
        z += n.z;
        return *this;
    }
    Normal3<T> operator-(const Normal3<T>& n) const {
        //DCHECK(!n.HasNaNs());
        return Normal3<T>(x - n.x, y - n.y, z - n.z);
    }

    Normal3<T>& operator-=(const Normal3<T>& n) {
        //DCHECK(!n.HasNaNs());
        x -= n.x;
        y -= n.y;
        z -= n.z;
        return *this;
    }
    bool HasNaNs() const { return isNaN(x) || isNaN(y) || isNaN(z); }
    template <typename U>
    Normal3<T> operator*(U f) const {
        return Normal3<T>(f * x, f * y, f * z);
    }

    template <typename U>
    Normal3<T>& operator*=(U f) {
        x *= f;
        y *= f;
        z *= f;
        return *this;
    }
    template <typename U>
    Normal3<T> operator/(U f) const {
        //CHECK_NE(f, 0);
        assert(f != 0);
        Float inv = (Float)1 / f;
        return Normal3<T>(x * inv, y * inv, z * inv);
    }

    template <typename U>
    Normal3<T>& operator/=(U f) {
        //CHECK_NE(f, 0);
        assert(f != 0);
        Float inv = (Float)1 / f;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }
    Float LengthSquared() const { return x * x + y * y + z * z; }
    Float Length() const { return std::sqrt(LengthSquared()); }

#ifndef NDEBUG
    Normal3<T>(const Normal3<T>& n) {
        DCHECK(!n.HasNaNs());
        x = n.x;
        y = n.y;
        z = n.z;
    }

    Normal3<T>& operator=(const Normal3<T>& n) {
        DCHECK(!n.HasNaNs());
        x = n.x;
        y = n.y;
        z = n.z;
        return *this;
    }
#endif  // !NDEBUG
    explicit Normal3<T>(const Vector3<T>& v) : x(v.x), y(v.y), z(v.z) {
        //DCHECK(!v.HasNaNs());
    }
    bool operator==(const Normal3<T>& n) const {
        return x == n.x && y == n.y && z == n.z;
    }
    bool operator!=(const Normal3<T>& n) const {
        return x != n.x || y != n.y || z != n.z;
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

    // Normal3 Public Data
    T x, y, z;
};

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Normal3<T>& v) {
    os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
    return os;
}

template <>
inline std::ostream& operator<<(std::ostream& os, const Normal3<Float>& v) {
    //os << StringPrintf("[ %f, %f, %f ]", v.x, v.y, v.z);
    os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
    return os;
}

//typedef Normal3<Float> Normal3f;
using Normal3f = Normal3<Float>;

//unity functions
template <typename T, typename U>
inline Vector3<T> operator*(U s, const Vector3<T>& v) {
    return v * s;
}
template <typename T>
Vector3<T> Abs(const Vector3<T>& v) {
    return Vector3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

template <typename T>
inline T Dot(const Vector3<T>& v1, const Vector3<T>& v2) {
    //DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <typename T>
inline T AbsDot(const Vector3<T>& v1, const Vector3<T>& v2) {
    //DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
    return std::abs(Dot(v1, v2));
}

template <typename T>
inline Vector3<T> Cross(const Vector3<T>& v1, const Vector3<T>& v2) {
    //DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
        (v1x * v2y) - (v1y * v2x));
}

template <typename T>
inline Vector3<T> Cross(const Vector3<T>& v1, const Normal3<T>& v2) {
    //DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
    double v1x = v1.x, v1y = v1.y, v1z = v1.z;
    double v2x = v2.x, v2y = v2.y, v2z = v2.z;
    return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
        (v1x * v2y) - (v1y * v2x));
}

template <typename T>
inline Vector3<T> Cross(const Normal3<T>& v1, const Vector3<T>& v2) {
    //DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
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
    if (std::abs(v1.x) > std::abs(v1.y))
        *v2 = Vector3<T>(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
    else
        *v2 = Vector3<T>(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);
    *v3 = Cross(v1, *v2);
}



template <typename T, typename U>
inline Vector2<T> operator*(U f, const Vector2<T>& v) {
    return v * f;
}
template <typename T>
inline Float Dot(const Vector2<T>& v1, const Vector2<T>& v2) {
    //DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
    return v1.x * v2.x + v1.y * v2.y;
}

template <typename T>
inline Float AbsDot(const Vector2<T>& v1, const Vector2<T>& v2) {
    //DCHECK(!v1.HasNaNs() && !v2.HasNaNs());
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

template <typename T>
inline Float Distance(const Point3<T>& p1, const Point3<T>& p2) {
    return (p1 - p2).Length();
}

template <typename T>
inline Float DistanceSquared(const Point3<T>& p1, const Point3<T>& p2) {
    return (p1 - p2).LengthSquared();
}

template <typename T, typename U>
inline Point3<T> operator*(U f, const Point3<T>& p) {
    //DCHECK(!p.HasNaNs());
    return p * f;
}

template <typename T>
Point3<T> Lerp(Float t, const Point3<T>& p0, const Point3<T>& p1) {
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
inline Float Distance(const Point2<T>& p1, const Point2<T>& p2) {
    return (p1 - p2).Length();
}

template <typename T>
inline Float DistanceSquared(const Point2<T>& p1, const Point2<T>& p2) {
    return (p1 - p2).LengthSquared();
}


template <typename T, typename U>
inline Point2<T> operator*(U f, const Point2<T>& p) {
    //DCHECK(!p.HasNaNs());
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
Point2<T> Lerp(Float t, const Point2<T>& v0, const Point2<T>& v1) {
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

template <typename T, typename U>
inline Normal3<T> operator*(U f, const Normal3<T>& n) {
    return Normal3<T>(f * n.x, f * n.y, f * n.z);
}

template <typename T>
inline Normal3<T> Normalize(const Normal3<T>& n) {
    return n / n.Length();
}


template <typename T>
inline T Dot(const Normal3<T>& n1, const Vector3<T>& v2) {
    //DCHECK(!n1.HasNaNs() && !v2.HasNaNs());
    return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
}

template <typename T>
inline T Dot(const Vector3<T>& v1, const Normal3<T>& n2) {
    //DCHECK(!v1.HasNaNs() && !n2.HasNaNs());
    return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}

template <typename T>
inline T Dot(const Normal3<T>& n1, const Normal3<T>& n2) {
    //DCHECK(!n1.HasNaNs() && !n2.HasNaNs());
    return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}

template <typename T>
inline T AbsDot(const Normal3<T>& n1, const Vector3<T>& v2) {
    //DCHECK(!n1.HasNaNs() && !v2.HasNaNs());
    return std::abs(n1.x * v2.x + n1.y * v2.y + n1.z * v2.z);
}

template <typename T>
inline T AbsDot(const Vector3<T>& v1, const Normal3<T>& n2) {
    //DCHECK(!v1.HasNaNs() && !n2.HasNaNs());
    return std::abs(v1.x * n2.x + v1.y * n2.y + v1.z * n2.z);
}

template <typename T>
inline T AbsDot(const Normal3<T>& n1, const Normal3<T>& n2) {
    //DCHECK(!n1.HasNaNs() && !n2.HasNaNs());
    return std::abs(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
}

template <typename T>
inline Normal3<T> Faceforward(const Normal3<T>& n, const Vector3<T>& v) {
    return (Dot(n, v) < 0.f) ? -n : n;
}

template <typename T>
inline Normal3<T> Faceforward(const Normal3<T>& n, const Normal3<T>& n2) {
    return (Dot(n, n2) < 0.f) ? -n : n;
}

template <typename T>
inline Vector3<T> Faceforward(const Vector3<T>& v, const Vector3<T>& v2) {
    return (Dot(v, v2) < 0.f) ? -v : v;
}

template <typename T>
inline Vector3<T> Faceforward(const Vector3<T>& v, const Normal3<T>& n2) {
    return (Dot(v, n2) < 0.f) ? -v : v;
}

template <typename T>
Normal3<T> Abs(const Normal3<T>& v) {
    return Normal3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

inline Point3f OffsetRayOrigin(const Point3f& p, const Vector3f& pError,
    const Normal3f& n, const Vector3f& w) {
    Float d = Dot(Abs(n), pError);
    Vector3f offset = d * Vector3f(n);
    if (Dot(w, n) < 0) offset = -offset;
    Point3f po = p + offset;
    // Round offset point _po_ away from _p_
    for (int i = 0; i < 3; ++i) {
        if (offset[i] > 0)
            po[i] = NextFloatUp(po[i]);
        else if (offset[i] < 0)
            po[i] = NextFloatDown(po[i]);
    }
    return po;
}

inline Vector3f SphericalDirection(Float sinTheta, Float cosTheta, Float phi) {
    return Vector3f(sinTheta * std::cos(phi), sinTheta * std::sin(phi),
        cosTheta);
}

inline Vector3f SphericalDirection(Float sinTheta, Float cosTheta, Float phi,
    const Vector3f& x, const Vector3f& y,
    const Vector3f& z) {
    return sinTheta * std::cos(phi) * x + sinTheta * std::sin(phi) * y +
        cosTheta * z;
}

inline Float SphericalTheta(const Vector3f& v) {
    return std::acos(Clamp(v.z, -1, 1));
}

inline Float SphericalPhi(const Vector3f& v) {
    Float p = std::atan2(v.y, v.x);
    return (p < 0) ? (p + 2 * pi) : p;
}


// vec3 Utility Functions
/*
inline std::ostream& operator<<(std::ostream& out, const vec3& v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline vec3 operator+(const vec3& u, const vec3& v) {
    return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vec3 operator-(const vec3& u, const vec3& v) {
    return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vec3 operator*(const vec3& u, const vec3& v) {
    return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vec3 operator*(double t, const vec3& v) {
    return vec3(t * v.e[0], t * v.e[1], t * v.e[2]);
}

inline vec3 operator*(const vec3& v, double t) {
    return t * v;
}

inline vec3 operator/(vec3 v, double t) {
    return (1 / t) * v;
}

inline double dot(const vec3& u, const vec3& v) {
    return u.e[0] * v.e[0]
        + u.e[1] * v.e[1]
        + u.e[2] * v.e[2];
}

inline vec3 cross(const vec3& u, const vec3& v) {
    return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
        u.e[2] * v.e[0] - u.e[0] * v.e[2],
        u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

*/

//origin
//template <typename T>
inline Vector3f unit_vector(Vector3f v) {
    return v / v.Length();
}

//template <typename T>
inline Vector3f random_in_unit_sphere() {
    while (true) {
        auto p = Vector3f::random(-1, 1);
        if (p.LengthSquared() >= 1) continue;
        return p;
    }
}

inline Vector3f random_unit_vector() {
    return unit_vector(random_in_unit_sphere());
}

inline Vector3f random_in_hemisphere(const Normal3f& normal) {
    Vector3f in_unit_sphere = random_in_unit_sphere();
    if (Dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

inline Vector3f random_in_unit_disk() {
    while (true) {
        auto p = Vector3f(random_Float(-1, 1), random_Float(-1, 1), 0);
        if (p.LengthSquared() >= 1) continue;
        return p;
    }
}

inline Vector3f reflect(const Vector3f& v, const Normal3f& n) {
    return v - 2 * Dot(v, n) * Vector3f(n);
}

inline Vector3f refract(const Vector3f& uv, const Normal3f& n, double etai_over_etat) {
    auto cos_theta = fmin(Dot(-uv, n), 1.0);
    Vector3f r_out_perp = etai_over_etat * (uv + cos_theta * Vector3f(n));
    Vector3f r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.LengthSquared())) * Vector3f(n);
    return r_out_perp + r_out_parallel;
}

inline Vector3f random_cosine_direction() {
    auto r1 = random_Float();
    auto r2 = random_Float();
    auto z = sqrt(1 - r2);

    auto phi = 2 * pi * r1;
    auto x = cos(phi) * sqrt(r2);
    auto y = sin(phi) * sqrt(r2);

    return Vector3f(x, y, z);
}

#endif