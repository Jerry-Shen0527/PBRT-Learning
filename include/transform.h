#pragma once

#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "rtweekend.h"
#include "aabb.h"
#include "hittable.h"

// Matrix4x4 Declarations
struct Matrix4x4 {
    // Matrix4x4 Public Methods
    Matrix4x4() {
        m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.f;
        m[0][1] = m[0][2] = m[0][3] = m[1][0] = m[1][2] = m[1][3] = m[2][0] =
            m[2][1] = m[2][3] = m[3][0] = m[3][1] = m[3][2] = 0.f;
    }
    Matrix4x4(Float mat[4][4]);
    Matrix4x4(Float t00, Float t01, Float t02, Float t03, Float t10, Float t11,
        Float t12, Float t13, Float t20, Float t21, Float t22, Float t23,
        Float t30, Float t31, Float t32, Float t33);
    bool operator==(const Matrix4x4& m2) const {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (m[i][j] != m2.m[i][j]) return false;
        return true;
    }
    bool operator!=(const Matrix4x4& m2) const {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (m[i][j] != m2.m[i][j]) return true;
        return false;
    }
    friend Matrix4x4 Transpose(const Matrix4x4&);
    void Print(FILE* f) const {
        fprintf(f, "[ ");
        for (int i = 0; i < 4; ++i) {
            fprintf(f, "  [ ");
            for (int j = 0; j < 4; ++j) {
                fprintf(f, "%f", m[i][j]);
                if (j != 3) fprintf(f, ", ");
            }
            fprintf(f, " ]\n");
        }
        fprintf(f, " ] ");
    }
    static Matrix4x4 Mul(const Matrix4x4& m1, const Matrix4x4& m2) {
        Matrix4x4 r;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                r.m[i][j] = m1.m[i][0] * m2.m[0][j] + m1.m[i][1] * m2.m[1][j] +
                m1.m[i][2] * m2.m[2][j] + m1.m[i][3] * m2.m[3][j];
        return r;
    }
    friend Matrix4x4 Inverse(const Matrix4x4&);

    friend std::ostream& operator<<(std::ostream& os, const Matrix4x4& m) {
        // clang-format off
        os << "[" <<
            m.m[0][0] << " " << m.m[0][1] << " " << m.m[0][2] << " " << m.m[0][3] << "]\n" <<
            m.m[1][0] << " " << m.m[1][1] << " " << m.m[1][2] << " " << m.m[1][3] << "]\n" <<
            m.m[2][0] << " " << m.m[2][1] << " " << m.m[2][2] << " " << m.m[2][3] << "]\n" <<
            m.m[3][0] << " " << m.m[3][1] << " " << m.m[3][2] << " " << m.m[3][3] << "]\n";
        // clang-format on
        return os;
    }

    Float m[4][4];
};
inline Matrix4x4::Matrix4x4(Float mat[4][4]) { memcpy(m, mat, 16 * sizeof(Float)); }

inline Matrix4x4::Matrix4x4(Float t00, Float t01, Float t02, Float t03, Float t10,
    Float t11, Float t12, Float t13, Float t20, Float t21,
    Float t22, Float t23, Float t30, Float t31, Float t32,
    Float t33) {
    m[0][0] = t00;
    m[0][1] = t01;
    m[0][2] = t02;
    m[0][3] = t03;
    m[1][0] = t10;
    m[1][1] = t11;
    m[1][2] = t12;
    m[1][3] = t13;
    m[2][0] = t20;
    m[2][1] = t21;
    m[2][2] = t22;
    m[2][3] = t23;
    m[3][0] = t30;
    m[3][1] = t31;
    m[3][2] = t32;
    m[3][3] = t33;
}

inline Matrix4x4 Transpose(const Matrix4x4& m) {
    return Matrix4x4(m.m[0][0], m.m[1][0], m.m[2][0], m.m[3][0], m.m[0][1],
        m.m[1][1], m.m[2][1], m.m[3][1], m.m[0][2], m.m[1][2],
        m.m[2][2], m.m[3][2], m.m[0][3], m.m[1][3], m.m[2][3],
        m.m[3][3]);
}

inline Matrix4x4 Inverse(const Matrix4x4& m) {
    int indxc[4], indxr[4];
    int ipiv[4] = { 0, 0, 0, 0 };
    Float minv[4][4];
    memcpy(minv, m.m, 4 * 4 * sizeof(Float));
    for (int i = 0; i < 4; i++) {
        int irow = 0, icol = 0;
        Float big = 0.f;
        // Choose pivot
        for (int j = 0; j < 4; j++) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < 4; k++) {
                    if (ipiv[k] == 0) {
                        if (std::abs(minv[j][k]) >= big) {
                            big = Float(std::abs(minv[j][k]));
                            irow = j;
                            icol = k;
                        }
                    }
                    else if (ipiv[k] > 1)
                        std::cerr<<"Singular matrix in MatrixInvert"<<std::endl;
                }
            }
        }
        ++ipiv[icol];
        // Swap rows _irow_ and _icol_ for pivot
        if (irow != icol) {
            for (int k = 0; k < 4; ++k) std::swap(minv[irow][k], minv[icol][k]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (minv[icol][icol] == 0.f) std::cerr << "Singular matrix in MatrixInvert" << std::endl;

        // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
        Float pivinv = 1. / minv[icol][icol];
        minv[icol][icol] = 1.;
        for (int j = 0; j < 4; j++) minv[icol][j] *= pivinv;

        // Subtract this row from others to zero out their columns
        for (int j = 0; j < 4; j++) {
            if (j != icol) {
                Float save = minv[j][icol];
                minv[j][icol] = 0;
                for (int k = 0; k < 4; k++) minv[j][k] -= minv[icol][k] * save;
            }
        }
    }
    // Swap columns to reflect permutation
    for (int j = 3; j >= 0; j--) {
        if (indxr[j] != indxc[j]) {
            for (int k = 0; k < 4; k++)
                std::swap(minv[k][indxr[j]], minv[k][indxc[j]]);
        }
    }
    return Matrix4x4(minv);
}

// Transform Declarations
class Transform {
public:
    // Transform Public Methods
    Transform() {}
    Transform(const Float mat[4][4]) {
        m = Matrix4x4(mat[0][0], mat[0][1], mat[0][2], mat[0][3], mat[1][0],
            mat[1][1], mat[1][2], mat[1][3], mat[2][0], mat[2][1],
            mat[2][2], mat[2][3], mat[3][0], mat[3][1], mat[3][2],
            mat[3][3]);
        mInv = Inverse(m);
    }
    Transform(const Matrix4x4& m) : m(m), mInv(Inverse(m)) {}
    Transform(const Matrix4x4& m, const Matrix4x4& mInv) : m(m), mInv(mInv) {}
    //void Print(FILE* f) const;
    friend Transform Inverse(const Transform& t) {
        return Transform(t.mInv, t.m);
    }
    friend Transform Transpose(const Transform& t) {
        return Transform(Transpose(t.m), Transpose(t.mInv));
    }
    bool operator==(const Transform& t) const {
        return t.m == m && t.mInv == mInv;
    }
    bool operator!=(const Transform& t) const {
        return t.m != m || t.mInv != mInv;
    }
    bool operator<(const Transform& t2) const {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) {
                if (m.m[i][j] < t2.m.m[i][j]) return true;
                if (m.m[i][j] > t2.m.m[i][j]) return false;
            }
        return false;
    }
    bool IsIdentity() const {
        return (m.m[0][0] == 1.f && m.m[0][1] == 0.f && m.m[0][2] == 0.f &&
            m.m[0][3] == 0.f && m.m[1][0] == 0.f && m.m[1][1] == 1.f &&
            m.m[1][2] == 0.f && m.m[1][3] == 0.f && m.m[2][0] == 0.f &&
            m.m[2][1] == 0.f && m.m[2][2] == 1.f && m.m[2][3] == 0.f &&
            m.m[3][0] == 0.f && m.m[3][1] == 0.f && m.m[3][2] == 0.f &&
            m.m[3][3] == 1.f);
    }
    const Matrix4x4& GetMatrix() const { return m; }
    const Matrix4x4& GetInverseMatrix() const { return mInv; }
    bool HasScale() const {
        Float la2 = (*this)(vec3(1, 0, 0)).LengthSquared();
        Float lb2 = (*this)(vec3(0, 1, 0)).LengthSquared();
        Float lc2 = (*this)(vec3(0, 0, 1)).LengthSquared();
#define NOT_ONE(x) ((x) < .999f || (x) > 1.001f)
        return (NOT_ONE(la2) || NOT_ONE(lb2) || NOT_ONE(lc2));
#undef NOT_ONE
    }
    inline point3 operator()(const point3& p) const;
    inline vec3 operator()(const vec3& v) const;
    inline ray operator()(const ray& r) const;
    inline Normal operator()(const Normal& n) const;
    //inline RayDifferential operator()(const RayDifferential& r) const;
    aabb operator()(const aabb& b) const;
    Transform operator*(const Transform& t2) const;
    bool SwapsHandedness() const;
    SurfaceInteraction operator()(const SurfaceInteraction& si) const;
    inline point3 operator()(const point3& pt,
        vec3* absError) const;
    inline point3 operator()(const point3& p, const vec3& pError,
        vec3* pTransError) const;
    inline vec3 operator()(const vec3& v,
        vec3* vTransError) const;
    inline vec3 operator()(const vec3& v, const vec3& vError,
        vec3* vTransError) const;
    inline ray operator()(const ray& r, vec3* oError,
        vec3* dError) const;
    inline ray operator()(const ray& r, const vec3& oErrorIn,
        const vec3& dErrorIn, vec3* oErrorOut,
        vec3* dErrorOut) const;

    friend std::ostream& operator<<(std::ostream& os, const Transform& t) {
        os << "t=" << t.m << ", inv=" << t.mInv;
        return os;
    }

private:
    // Transform Private Data
    Matrix4x4 m, mInv;
    friend class AnimatedTransform;
    friend struct Quaternion;
};


Transform Translate(const vec3& delta);
Transform Scale(Float x, Float y, Float z);
Transform RotateX(Float theta);
Transform RotateY(Float theta);
Transform RotateZ(Float theta);
Transform Rotate(Float theta, const vec3& axis);
Transform LookAt(const point3& pos, const point3& look, const vec3& up);
Transform Orthographic(Float znear, Float zfar);
Transform Perspective(Float fov, Float znear, Float zfar);
bool SolveLinearSystem2x2(const Float A[2][2], const Float B[2], Float* x0,
    Float* x1);

// Transform Inline Functions
inline point3 Transform::operator()(const point3& p) const {
    Float x = p.x, y = p.y, z = p.z;
    Float xp = m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z + m.m[0][3];
    Float yp = m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z + m.m[1][3];
    Float zp = m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z + m.m[2][3];
    Float wp = m.m[3][0] * x + m.m[3][1] * y + m.m[3][2] * z + m.m[3][3];
    ZERO_DENOMINATOR(wp);
    if (wp == 1)
        return point3(xp, yp, zp);
    else
        return point3(xp, yp, zp) / wp;
}

inline vec3 Transform::operator()(const vec3& v) const {
    Float x = v.x, y = v.y, z = v.z;
    return vec3(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
        m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
        m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

inline Normal Transform::operator()(const Normal & n) const {
    Float x = n.x, y = n.y, z = n.z;
    return Normal(mInv.m[0][0] * x + mInv.m[1][0] * y + mInv.m[2][0] * z,
        mInv.m[0][1] * x + mInv.m[1][1] * y + mInv.m[2][1] * z,
        mInv.m[0][2] * x + mInv.m[1][2] * y + mInv.m[2][2] * z);
}

inline ray Transform::operator()(const ray& r) const {
    vec3 oError;
    point3 o = (*this)(r.o, &oError);
    vec3 d = (*this)(r.d);
    // Offset ray origin to edge of error bounds and compute _tMax_
    Float lengthSquared = d.LengthSquared();
    Float tMax = r.tMax;
    if (lengthSquared > 0) {
        Float dt = Dot(Abs(d), oError) / lengthSquared;
        o += d * dt;
        tMax -= dt;
    }
    //return ray(o, d, tMax, r.time, r.medium);
    return ray(o, d, r.time, r.tMin, r.tMax, r.depth);
}
/*
inline RayDifferential Transform::operator()(const RayDifferential& r) const {
    Ray tr = (*this)(ray(r));
    RayDifferential ret(tr.o, tr.d, tr.tMax, tr.time, tr.medium);
    ret.hasDifferentials = r.hasDifferentials;
    ret.rxOrigin = (*this)(r.rxOrigin);
    ret.ryOrigin = (*this)(r.ryOrigin);
    ret.rxDirection = (*this)(r.rxDirection);
    ret.ryDirection = (*this)(r.ryDirection);
    return ret;
}*/

inline point3 Transform::operator()(const point3& p,
    vec3* pError) const {
    Float x = p.x, y = p.y, z = p.z;
    // Compute transformed coordinates from point _pt_
    Float xp = (m.m[0][0] * x + m.m[0][1] * y) + (m.m[0][2] * z + m.m[0][3]);
    Float yp = (m.m[1][0] * x + m.m[1][1] * y) + (m.m[1][2] * z + m.m[1][3]);
    Float zp = (m.m[2][0] * x + m.m[2][1] * y) + (m.m[2][2] * z + m.m[2][3]);
    Float wp = (m.m[3][0] * x + m.m[3][1] * y) + (m.m[3][2] * z + m.m[3][3]);

    // Compute absolute error for transformed point
    Float xAbsSum = (std::abs(m.m[0][0] * x) + std::abs(m.m[0][1] * y) +
        std::abs(m.m[0][2] * z) + std::abs(m.m[0][3]));
    Float yAbsSum = (std::abs(m.m[1][0] * x) + std::abs(m.m[1][1] * y) +
        std::abs(m.m[1][2] * z) + std::abs(m.m[1][3]));
    Float zAbsSum = (std::abs(m.m[2][0] * x) + std::abs(m.m[2][1] * y) +
        std::abs(m.m[2][2] * z) + std::abs(m.m[2][3]));
    *pError = gamma(3) * vec3(xAbsSum, yAbsSum, zAbsSum);
    ZERO_DENOMINATOR(wp);
    if (wp == 1)
        return point3(xp, yp, zp);
    else
        return point3(xp, yp, zp) / wp;
}

inline point3 Transform::operator()(const point3& pt,
    const vec3& ptError,
    vec3* absError) const {
    Float x = pt.x, y = pt.y, z = pt.z;
    Float xp = (m.m[0][0] * x + m.m[0][1] * y) + (m.m[0][2] * z + m.m[0][3]);
    Float yp = (m.m[1][0] * x + m.m[1][1] * y) + (m.m[1][2] * z + m.m[1][3]);
    Float zp = (m.m[2][0] * x + m.m[2][1] * y) + (m.m[2][2] * z + m.m[2][3]);
    Float wp = (m.m[3][0] * x + m.m[3][1] * y) + (m.m[3][2] * z + m.m[3][3]);
    absError->x =
        (gamma(3) + (Float)1) *
        (std::abs(m.m[0][0]) * ptError.x + std::abs(m.m[0][1]) * ptError.y +
            std::abs(m.m[0][2]) * ptError.z) +
        gamma(3) * (std::abs(m.m[0][0] * x) + std::abs(m.m[0][1] * y) +
            std::abs(m.m[0][2] * z) + std::abs(m.m[0][3]));
    absError->y =
        (gamma(3) + (Float)1) *
        (std::abs(m.m[1][0]) * ptError.x + std::abs(m.m[1][1]) * ptError.y +
            std::abs(m.m[1][2]) * ptError.z) +
        gamma(3) * (std::abs(m.m[1][0] * x) + std::abs(m.m[1][1] * y) +
            std::abs(m.m[1][2] * z) + std::abs(m.m[1][3]));
    absError->z =
        (gamma(3) + (Float)1) *
        (std::abs(m.m[2][0]) * ptError.x + std::abs(m.m[2][1]) * ptError.y +
            std::abs(m.m[2][2]) * ptError.z) +
        gamma(3) * (std::abs(m.m[2][0] * x) + std::abs(m.m[2][1] * y) +
            std::abs(m.m[2][2] * z) + std::abs(m.m[2][3]));
    ZERO_DENOMINATOR(wp);
    if (wp == 1.)
        return point3(xp, yp, zp);
    else
        return point3(xp, yp, zp) / wp;
}

inline vec3 Transform::operator()(const vec3& v,
    vec3* absError) const {
    Float x = v.x, y = v.y, z = v.z;
    absError->x =
        gamma(3) * (std::abs(m.m[0][0] * v.x) + std::abs(m.m[0][1] * v.y) +
            std::abs(m.m[0][2] * v.z));
    absError->y =
        gamma(3) * (std::abs(m.m[1][0] * v.x) + std::abs(m.m[1][1] * v.y) +
            std::abs(m.m[1][2] * v.z));
    absError->z =
        gamma(3) * (std::abs(m.m[2][0] * v.x) + std::abs(m.m[2][1] * v.y) +
            std::abs(m.m[2][2] * v.z));
    return vec3(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
        m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
        m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

inline vec3 Transform::operator()(const vec3& v,
    const vec3& vError,
    vec3* absError) const {
    Float x = v.x, y = v.y, z = v.z;
    absError->x =
        (gamma(3) + (Float)1) *
        (std::abs(m.m[0][0]) * vError.x + std::abs(m.m[0][1]) * vError.y +
            std::abs(m.m[0][2]) * vError.z) +
        gamma(3) * (std::abs(m.m[0][0] * v.x) + std::abs(m.m[0][1] * v.y) +
            std::abs(m.m[0][2] * v.z));
    absError->y =
        (gamma(3) + (Float)1) *
        (std::abs(m.m[1][0]) * vError.x + std::abs(m.m[1][1]) * vError.y +
            std::abs(m.m[1][2]) * vError.z) +
        gamma(3) * (std::abs(m.m[1][0] * v.x) + std::abs(m.m[1][1] * v.y) +
            std::abs(m.m[1][2] * v.z));
    absError->z =
        (gamma(3) + (Float)1) *
        (std::abs(m.m[2][0]) * vError.x + std::abs(m.m[2][1]) * vError.y +
            std::abs(m.m[2][2]) * vError.z) +
        gamma(3) * (std::abs(m.m[2][0] * v.x) + std::abs(m.m[2][1] * v.y) +
            std::abs(m.m[2][2] * v.z));
    return vec3(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
        m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
        m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

inline ray Transform::operator()(const ray& r, vec3* oError,
    vec3* dError) const {
    point3 o = (*this)(r.o, oError);
    vec3 d = (*this)(r.d, dError);
    Float tMax = r.tMax;
    Float lengthSquared = d.LengthSquared();
    if (lengthSquared > 0) {
        Float dt = Dot(Abs(d), *oError) / lengthSquared;
        o += d * dt;
        //        tMax -= dt;
    }
    return ray(o, d, r.time, r.tMin, r.tMax, r.depth);
}

inline ray Transform::operator()(const ray& r, const vec3& oErrorIn,
    const vec3& dErrorIn, vec3* oErrorOut,
    vec3* dErrorOut) const {
    point3 o = (*this)(r.o, oErrorIn, oErrorOut);
    vec3 d = (*this)(r.d, dErrorIn, dErrorOut);
    Float tMax = r.tMax;
    Float lengthSquared = d.LengthSquared();
    if (lengthSquared > 0) {
        Float dt = Dot(Abs(d), *oErrorOut) / lengthSquared;
        o += d * dt;
        //        tMax -= dt;
    }
    return ray(o, d, r.time, r.tMin, r.tMax, r.depth);
}



#endif // !TRANSFORM_H
