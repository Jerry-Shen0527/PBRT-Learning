#ifndef MTRANSFORM_H
#define MTRANSFORM_H
#include <iostream>
#include "vec3.h"
#include "ray.h"
#include <assert.h>
#include "rtweekend.h"

typedef double Float;
struct Matrix4x4 {
    // Matrix4x4 Public Methods
    Matrix4x4() {
        m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.0;
        m[0][1] = m[0][2] = m[0][3] = m[1][0] = m[1][2] = m[1][3] = m[2][0] =
            m[2][1] = m[2][3] = m[3][0] = m[3][1] = m[3][2] = 0.0;
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

    //friend std::ostream& operator<<(std::ostream& os, const Matrix4x4& m) {
    //    // clang-format off
    //    os << StringPrintf("[ [ %f, %f, %f, %f ] "
    //        "[ %f, %f, %f, %f ] "
    //        "[ %f, %f, %f, %f ] "
    //        "[ %f, %f, %f, %f ] ]",
    //        m.m[0][0], m.m[0][1], m.m[0][2], m.m[0][3],
    //        m.m[1][0], m.m[1][1], m.m[1][2], m.m[1][3],
    //        m.m[2][0], m.m[2][1], m.m[2][2], m.m[2][3],
    //        m.m[3][0], m.m[3][1], m.m[3][2], m.m[3][3]);
    //    // clang-format on
    //    return os;
    //}

    Float m[4][4];
};

// Transform Declarations
class Transform {
public:
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

    inline vec3 vector_transform(const vec3& v) const;
    inline point3 point_transform(const point3& v) const;
    inline vec3 times_inverse(const vec3& v) const;

    inline ray ray_transform(const ray& r) const;


    Transform operator*(const Transform& t2) const;
    



private:
    // Transform Private Data
    Matrix4x4 m, mInv;
    
};





inline vec3 Transform::vector_transform(const vec3& v) const {
    double x = v.e[0], y = v.e[1], z = v.e[2];
    return vec3(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
        m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
        m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

inline point3 Transform::point_transform(const point3& v) const
{
    double x = v.e[0], y = v.e[1], z = v.e[2];
    double xp = m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z + m.m[0][3];
    double yp = m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z + m.m[1][3];
    double zp = m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z + m.m[2][3];
    double wp = m.m[3][0] * x + m.m[3][1] * y + m.m[3][2] * z + m.m[3][3];
    assert(wp != 0);
    if (wp == 1)
        return point3(xp, yp, zp);
    else
        return point3(xp, yp, zp) / wp;
}

inline vec3 Transform::times_inverse(const vec3& v) const
{
    double x = v.e[0], y = v.e[1], z = v.e[2];
    return vec3(mInv.m[0][0] * x + mInv.m[1][0] * y + mInv.m[2][0] * z,
        mInv.m[0][1] * x + mInv.m[1][1] * y + mInv.m[2][1] * z,
        mInv.m[0][2] * x + mInv.m[1][2] * y + mInv.m[2][2] * z);
}

inline ray Transform::ray_transform(const ray& r) const
{
    ray t_ray;
    t_ray.orig=point_transform(r.orig);
    t_ray.dir = vector_transform(r.dir);
    t_ray.tm = r.tm;
    return t_ray;
}



Matrix4x4::Matrix4x4(Float mat[4][4]) { memcpy(m, mat, 16 * sizeof(Float)); }

Matrix4x4::Matrix4x4(Float t00, Float t01, Float t02, Float t03, Float t10,
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

Matrix4x4 Transpose(const Matrix4x4& m) {
    return Matrix4x4(m.m[0][0], m.m[1][0], m.m[2][0], m.m[3][0], m.m[0][1],
        m.m[1][1], m.m[2][1], m.m[3][1], m.m[0][2], m.m[1][2],
        m.m[2][2], m.m[3][2], m.m[0][3], m.m[1][3], m.m[2][3],
        m.m[3][3]);
}

inline void Error(const char* a)
{
    std::cerr << a;
}

Matrix4x4 Inverse(const Matrix4x4& m) {
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
                        Error("Singular matrix in MatrixInvert");
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
        if (minv[icol][icol] == 0.f) Error("Singular matrix in MatrixInvert");

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

Transform Translate(const vec3& delta);
Transform RotateX(Float theta);
Transform RotateY(Float theta);
Transform RotateZ(Float theta);
Transform Rotate(Float theta, const vec3& axis);





Transform Translate(const vec3& delta) {
    Matrix4x4 m(1, 0, 0, delta.e[0], 0, 1, 0, delta.e[1], 0, 0, 1, delta.e[2], 0, 0, 0,
        1);
    Matrix4x4 minv(1, 0, 0, -delta.e[0], 0, 1, 0, -delta.e[1], 0, 0, 1, -delta.e[2], 0,
        0, 0, 1);
    return Transform(m, minv);
}

Transform RotateX(Float theta) {
    Float sinTheta = std::sin(degrees_to_radians(theta));
    Float cosTheta = std::cos(degrees_to_radians(theta));
    Matrix4x4 m(1, 0, 0, 0, 0, cosTheta, -sinTheta, 0, 0, sinTheta, cosTheta, 0,
        0, 0, 0, 1);
    return Transform(m, Transpose(m));
}

Transform RotateY(Float theta) {
    Float sinTheta = std::sin(degrees_to_radians(theta));
    Float cosTheta = std::cos(degrees_to_radians(theta));
    Matrix4x4 m(cosTheta, 0, sinTheta, 0, 0, 1, 0, 0, -sinTheta, 0, cosTheta, 0,
        0, 0, 0, 1);
    return Transform(m, Transpose(m));
}

Transform RotateZ(Float theta) {
    Float sinTheta = std::sin(degrees_to_radians(theta));
    Float cosTheta = std::cos(degrees_to_radians(theta));
    Matrix4x4 m(cosTheta, -sinTheta, 0, 0, sinTheta, cosTheta, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 1);
    return Transform(m, Transpose(m));
}

Transform Rotate(Float theta, const vec3& axis) {
    vec3 a = unit_vector(axis);
    Float sinTheta = std::sin(degrees_to_radians(theta));
    Float cosTheta = std::cos(degrees_to_radians(theta));
    Matrix4x4 m;
    // Compute rotation of first basis vector
    m.m[0][0] = a.e[0] * a.e[0] + (1 - a.e[0] * a.e[0]) * cosTheta;
    m.m[0][1] = a.e[0] * a.e[1] * (1 - cosTheta) - a.e[2] * sinTheta;
    m.m[0][2] = a.e[0] * a.e[2] * (1 - cosTheta) + a.e[1] * sinTheta;
    m.m[0][3] = 0;

    // Compute rotations of second and third basis vectors
    m.m[1][0] = a.e[0] * a.e[1] * (1 - cosTheta) + a.e[2] * sinTheta;
    m.m[1][1] = a.e[1] * a.e[1] + (1 - a.e[1] * a.e[1]) * cosTheta;
    m.m[1][2] = a.e[1] * a.e[2] * (1 - cosTheta) - a.e[0] * sinTheta;
    m.m[1][3] = 0;

    m.m[2][0] = a.e[0] * a.e[2] * (1 - cosTheta) - a.e[1] * sinTheta;
    m.m[2][1] = a.e[1] * a.e[2] * (1 - cosTheta) + a.e[0] * sinTheta;
    m.m[2][2] = a.e[2] * a.e[2] + (1 - a.e[2] * a.e[2]) * cosTheta;
    m.m[2][3] = 0;
    return Transform(m, Transpose(m));
}

Transform Transform::operator*(const Transform& t2) const {
    return Transform(Matrix4x4::Mul(m, t2.m), Matrix4x4::Mul(t2.mInv, mInv));
}
#endif
