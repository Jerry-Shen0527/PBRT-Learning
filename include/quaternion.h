#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef QUATERNION_H
#define QUATERNION_H

#include "vec3.h"

class Transform;
// Quaternion Declarations
struct Quaternion {
    // Quaternion Public Methods
    Quaternion() : v(0, 0, 0), w(1) {}
    Quaternion& operator+=(const Quaternion& q) {
        v += q.v;
        w += q.w;
        return *this;
    }
    friend Quaternion operator+(const Quaternion& q1, const Quaternion& q2) {
        Quaternion ret = q1;
        return ret += q2;
    }
    Quaternion& operator-=(const Quaternion& q) {
        v -= q.v;
        w -= q.w;
        return *this;
    }
    Quaternion operator-() const {
        Quaternion ret;
        ret.v = -v;
        ret.w = -w;
        return ret;
    }
    friend Quaternion operator-(const Quaternion& q1, const Quaternion& q2) {
        Quaternion ret = q1;
        return ret -= q2;
    }
    Quaternion& operator*=(Float f) {
        v *= f;
        w *= f;
        return *this;
    }
    Quaternion operator*(Float f) const {
        Quaternion ret = *this;
        ret.v *= f;
        ret.w *= f;
        return ret;
    }
    Quaternion& operator/=(Float f) {
        v /= f;
        w /= f;
        return *this;
    }
    Quaternion operator/(Float f) const {
        Quaternion ret = *this;
        ret.v /= f;
        ret.w /= f;
        return ret;
    }
    Transform ToTransform() const;
    Quaternion(const Transform& t);

    //friend std::ostream& operator<<(std::ostream& os, const Quaternion& q) {
    //    os << StringPrintf("[ %f, %f, %f, %f ]", q.v.x, q.v.y, q.v.z,
    //        q.w);
    //    return os;
    //}

    // Quaternion Public Data
    Vector3f v;
    Float w;
};

Quaternion Slerp(Float t, const Quaternion& q1, const Quaternion& q2);

// Quaternion Inline Functions
inline Quaternion operator*(Float f, const Quaternion& q) { return q * f; }

inline Float Dot(const Quaternion& q1, const Quaternion& q2) {
    return Dot(q1.v, q2.v) + q1.w * q2.w;
}

inline Quaternion Normalize(const Quaternion& q) {
    return q / std::sqrt(Dot(q, q));
}


#endif // ! QUATERNION_H
