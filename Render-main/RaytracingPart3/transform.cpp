#include "transform.h"
bool SolveLinearSystem2x2(const double A[2][2], const double B[2], double* x0,
    double* x1) {
    double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    if (std::abs(det) < 1e-10f) return false;
    *x0 = (A[1][1] * B[0] - A[0][1] * B[1]) / det;
    *x1 = (A[0][0] * B[1] - A[1][0] * B[0]) / det;
    if (std::isnan(*x0) || std::isnan(*x1)) return false;
    return true;
}

Matrix4x4::Matrix4x4(double mat[4][4]) { memcpy(m, mat, 16 * sizeof(double)); }

Matrix4x4::Matrix4x4(double t00, double t01, double t02, double t03, double t10,
    double t11, double t12, double t13, double t20, double t21,
    double t22, double t23, double t30, double t31, double t32,
    double t33) {
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

Matrix4x4 Inverse(const Matrix4x4& m) {
    int indxc[4], indxr[4];
    int ipiv[4] = { 0, 0, 0, 0 };
    double minv[4][4];
    memcpy(minv, m.m, 4 * 4 * sizeof(double));
    for (int i = 0; i < 4; i++) {
        int irow = 0, icol = 0;
        double big = 0.f;
        // Choose pivot
        for (int j = 0; j < 4; j++) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < 4; k++) {
                    if (ipiv[k] == 0) {
                        if (std::abs(minv[j][k]) >= big) {
                            big = double(std::abs(minv[j][k]));
                            irow = j;
                            icol = k;
                        }
                    }
                    else if (ipiv[k] > 1) {
                        perror("Singular matrix in MatrixInvert");
                        exit(EXIT_FAILURE);
                    }
                        
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
        if (minv[icol][icol] == 0.f) {
            perror("Singular matrix in MatrixInvert");
            exit(EXIT_FAILURE);
        }

        // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
        double pivinv = 1. / minv[icol][icol];
        minv[icol][icol] = 1.;
        for (int j = 0; j < 4; j++) minv[icol][j] *= pivinv;

        // Subtract this row from others to zero out their columns
        for (int j = 0; j < 4; j++) {
            if (j != icol) {
                double save = minv[j][icol];
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

ray ConvertRayTrans(const Transform& trans, const ray& r, Vector3f* oError,
    Vector3f* dError, double& tmax, double& tmin) {
    Point3f o = ConvertPTrans(trans,r.origin(), *oError);
    //cout << o.x << endl;
    Vector3f d = trans(r.direction(), dError);
    double tMax = tmax;
    double lengthSquared = d.LengthSquared();
    if (lengthSquared > 0) {
        double dt = Dot(Abs(d), *oError) / lengthSquared;
        o += d * dt;
        //        tMax -= dt;
    }
    return ray(o, d, r.time());
}

Point3f ConvertPTrans(const Transform& trans, const Point3f& p, Vector3f& pError) {
    Matrix4x4 m = trans.GetMatrix();
    double x = p.x, y = p.y, z = p.z;
    //cout << x << endl;
    // Compute transformed coordinates from point _pt_
    double xp = (m.m[0][0] * x + m.m[0][1] * y) + (m.m[0][2] * z + m.m[0][3]);
    double yp = (m.m[1][0] * x + m.m[1][1] * y) + (m.m[1][2] * z + m.m[1][3]);
    double zp = (m.m[2][0] * x + m.m[2][1] * y) + (m.m[2][2] * z + m.m[2][3]);
    double wp = (m.m[3][0] * x + m.m[3][1] * y) + (m.m[3][2] * z + m.m[3][3]);
    //cout << xp << " " << m.m[0][3] << endl;

    // Compute absolute error for transformed point
    double xAbsSum = (std::abs(m.m[0][0] * x) + std::abs(m.m[0][1] * y) +
        std::abs(m.m[0][2] * z) + std::abs(m.m[0][3]));
    double yAbsSum = (std::abs(m.m[1][0] * x) + std::abs(m.m[1][1] * y) +
        std::abs(m.m[1][2] * z) + std::abs(m.m[1][3]));
    double zAbsSum = (std::abs(m.m[2][0] * x) + std::abs(m.m[2][1] * y) +
        std::abs(m.m[2][2] * z) + std::abs(m.m[2][3]));
    pError = gamma(3) * Vector3f(xAbsSum, yAbsSum, zAbsSum);
    assert(wp != 0);
    if (wp == 1)
        return Point3f(xp, yp, zp);
    else
        return Point3f(xp, yp, zp) / wp;
}

Point3f ConvertPTrans(const Transform& trans, const Point3f& p) {
    Matrix4x4 m = trans.GetMatrix();
    double x = p.x, y = p.y, z = p.z;
    // Compute transformed coordinates from point _pt_
    double xp = (m.m[0][0] * x + m.m[0][1] * y) + (m.m[0][2] * z + m.m[0][3]);
    double yp = (m.m[1][0] * x + m.m[1][1] * y) + (m.m[1][2] * z + m.m[1][3]);
    double zp = (m.m[2][0] * x + m.m[2][1] * y) + (m.m[2][2] * z + m.m[2][3]);
    double wp = (m.m[3][0] * x + m.m[3][1] * y) + (m.m[3][2] * z + m.m[3][3]);

    
        return Point3f(xp, yp, zp);
    
}

Vector3f ConvertVTrans(const Transform& trans, const Vector3f& v, Vector3f& vError) {
    Matrix4x4 m = trans.GetMatrix();
    double x = v.x, y = v.y, z = v.z;
    vError.x =
        gamma(3) * (std::abs(m.m[0][0] * v.x) + std::abs(m.m[0][1] * v.y) +
            std::abs(m.m[0][2] * v.z));
    vError.y =
        gamma(3) * (std::abs(m.m[1][0] * v.x) + std::abs(m.m[1][1] * v.y) +
            std::abs(m.m[1][2] * v.z));
    vError.z =
        gamma(3) * (std::abs(m.m[2][0] * v.x) + std::abs(m.m[2][1] * v.y) +
            std::abs(m.m[2][2] * v.z));
    return Vector3f(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
        m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
        m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

Vector3f ConvertVTrans(const Transform& trans, const Vector3f& v) {
    Matrix4x4 mInv = trans.GetMatrix();
    double x = v.x, y = v.y, z = v.z;
    return Vector3f(mInv.m[0][0] * x + mInv.m[1][0] * y + mInv.m[2][0] * z,
        mInv.m[0][1] * x + mInv.m[1][1] * y + mInv.m[2][1] * z,
        mInv.m[0][2] * x + mInv.m[1][2] * y + mInv.m[2][2] * z);
}

void Transform::Print(FILE* f) const { m.Print(f); }

Transform Translate(const Vector3f& delta) {
    Matrix4x4 m(1, 0, 0, delta.x, 0, 1, 0, delta.y, 0, 0, 1, delta.z, 0, 0, 0,
        1);
    Matrix4x4 minv(1, 0, 0, -delta.x, 0, 1, 0, -delta.y, 0, 0, 1, -delta.z, 0,
        0, 0, 1);
    return Transform(m, minv);
}

Transform Scale(double x, double y, double z) {
    Matrix4x4 m(x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1);
    Matrix4x4 minv(1 / x, 0, 0, 0, 0, 1 / y, 0, 0, 0, 0, 1 / z, 0, 0, 0, 0, 1);
    return Transform(m, minv);
}

Transform RotateX(double theta) {
    double sinTheta = std::sin(Radians(theta));
    double cosTheta = std::cos(Radians(theta));
    Matrix4x4 m(1, 0, 0, 0, 0, cosTheta, -sinTheta, 0, 0, sinTheta, cosTheta, 0,
        0, 0, 0, 1);
    return Transform(m, Transpose(m));
}

Transform RotateY(double theta) {
    double sinTheta = std::sin(Radians(theta));
    double cosTheta = std::cos(Radians(theta));
    Matrix4x4 m(cosTheta, 0, sinTheta, 0, 0, 1, 0, 0, -sinTheta, 0, cosTheta, 0,
        0, 0, 0, 1);
    return Transform(m, Transpose(m));
}

Transform RotateZ(double theta) {
    double sinTheta = std::sin(Radians(theta));
    double cosTheta = std::cos(Radians(theta));
    Matrix4x4 m(cosTheta, -sinTheta, 0, 0, sinTheta, cosTheta, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 1);
    return Transform(m, Transpose(m));
}

Transform Rotate(double theta, const Vector3f& axis) {
    Vector3f a = Normalize(axis);
    double sinTheta = std::sin(Radians(theta));
    double cosTheta = std::cos(Radians(theta));
    Matrix4x4 m;
    // Compute rotation of first basis vector
    m.m[0][0] = a.x * a.x + (1 - a.x * a.x) * cosTheta;
    m.m[0][1] = a.x * a.y * (1 - cosTheta) - a.z * sinTheta;
    m.m[0][2] = a.x * a.z * (1 - cosTheta) + a.y * sinTheta;
    m.m[0][3] = 0;

    // Compute rotations of second and third basis vectors
    m.m[1][0] = a.x * a.y * (1 - cosTheta) + a.z * sinTheta;
    m.m[1][1] = a.y * a.y + (1 - a.y * a.y) * cosTheta;
    m.m[1][2] = a.y * a.z * (1 - cosTheta) - a.x * sinTheta;
    m.m[1][3] = 0;

    m.m[2][0] = a.x * a.z * (1 - cosTheta) - a.y * sinTheta;
    m.m[2][1] = a.y * a.z * (1 - cosTheta) + a.x * sinTheta;
    m.m[2][2] = a.z * a.z + (1 - a.z * a.z) * cosTheta;
    m.m[2][3] = 0;
    return Transform(m, Transpose(m));
}

Transform LookAt(const Point3f& pos, const Point3f& look, const Vector3f& up) {
    Matrix4x4 cameraToWorld;
    // Initialize fourth column of viewing matrix
    cameraToWorld.m[0][3] = pos.x;
    cameraToWorld.m[1][3] = pos.y;
    cameraToWorld.m[2][3] = pos.z;
    cameraToWorld.m[3][3] = 1;

    // Initialize first three columns of viewing matrix
    Vector3f dir = Normalize(look - pos);
    if (Cross(Normalize(up), dir).Length() == 0) {
        printf("\"up\" vector (%f, %f, %f) and viewing direction (%f, %f, %f) "
            "passed to LookAt are pointing in the same direction.  Using "
            "the identity transformation.",
            up.x, up.y, up.z, dir.x, dir.y, dir.z);
        exit(EXIT_FAILURE);
        
        return Transform();
    }
    Vector3f right = Normalize(Cross(Normalize(up), dir));
    Vector3f newUp = Cross(dir, right);
    cameraToWorld.m[0][0] = right.x;
    cameraToWorld.m[1][0] = right.y;
    cameraToWorld.m[2][0] = right.z;
    cameraToWorld.m[3][0] = 0.;
    cameraToWorld.m[0][1] = newUp.x;
    cameraToWorld.m[1][1] = newUp.y;
    cameraToWorld.m[2][1] = newUp.z;
    cameraToWorld.m[3][1] = 0.;
    cameraToWorld.m[0][2] = dir.x;
    cameraToWorld.m[1][2] = dir.y;
    cameraToWorld.m[2][2] = dir.z;
    cameraToWorld.m[3][2] = 0.;
    return Transform(Inverse(cameraToWorld), cameraToWorld);
}

/*Bounds3f Transform::operator()(const Bounds3f& b) const {
    const Transform& M = *this;
    Bounds3f ret(M(Point3f(b.pMin.x, b.pMin.y, b.pMin.z)));
    ret = Union(ret, M(Point3f(b.pMax.x, b.pMin.y, b.pMin.z)));
    ret = Union(ret, M(Point3f(b.pMin.x, b.pMax.y, b.pMin.z)));
    ret = Union(ret, M(Point3f(b.pMin.x, b.pMin.y, b.pMax.z)));
    ret = Union(ret, M(Point3f(b.pMin.x, b.pMax.y, b.pMax.z)));
    ret = Union(ret, M(Point3f(b.pMax.x, b.pMax.y, b.pMin.z)));
    ret = Union(ret, M(Point3f(b.pMax.x, b.pMin.y, b.pMax.z)));
    ret = Union(ret, M(Point3f(b.pMax.x, b.pMax.y, b.pMax.z)));
    return ret;
}*/


bool Transform::SwapsHandedness() const {
    double det = m.m[0][0] * (m.m[1][1] * m.m[2][2] - m.m[1][2] * m.m[2][1]) -
        m.m[0][1] * (m.m[1][0] * m.m[2][2] - m.m[1][2] * m.m[2][0]) +
        m.m[0][2] * (m.m[1][0] * m.m[2][1] - m.m[1][1] * m.m[2][0]);
    return det < 0;
}

/*AnimatedTransform::AnimatedTransform(const Transform* startTransform,
    double startTime,
    const Transform* endTransform,
    double endTime)
    : startTransform(startTransform),
    endTransform(endTransform),
    startTime(startTime),
    endTime(endTime),
    actuallyAnimated(*startTransform != *endTransform) {
    if (!actuallyAnimated)
        return;
    Decompose(startTransform->m, &T[0], &R[0], &S[0]);
    Decompose(endTransform->m, &T[1], &R[1], &S[1]);
    // Flip _R[1]_ if needed to select shortest path
    if (Dot(R[0], R[1]) < 0) R[1] = -R[1];
    hasRotation = Dot(R[0], R[1]) < 0.9995f;
    // Compute terms of motion derivative function
    if (hasRotation) {
        double cosTheta = Dot(R[0], R[1]);
        double theta = std::acos(check::Clamp(cosTheta, -1, 1));
        Quaternion qperp = Normalize(R[1] - R[0] * cosTheta);

        double t0x = T[0].x;
        double t0y = T[0].y;
        double t0z = T[0].z;
        double t1x = T[1].x;
        double t1y = T[1].y;
        double t1z = T[1].z;
        double q0x = R[0].v.x;
        double q0y = R[0].v.y;
        double q0z = R[0].v.z;
        double q0w = R[0].w;
        double qperpx = qperp.v.x;
        double qperpy = qperp.v.y;
        double qperpz = qperp.v.z;
        double qperpw = qperp.w;
        double s000 = S[0].m[0][0];
        double s001 = S[0].m[0][1];
        double s002 = S[0].m[0][2];
        double s010 = S[0].m[1][0];
        double s011 = S[0].m[1][1];
        double s012 = S[0].m[1][2];
        double s020 = S[0].m[2][0];
        double s021 = S[0].m[2][1];
        double s022 = S[0].m[2][2];
        double s100 = S[1].m[0][0];
        double s101 = S[1].m[0][1];
        double s102 = S[1].m[0][2];
        double s110 = S[1].m[1][0];
        double s111 = S[1].m[1][1];
        double s112 = S[1].m[1][2];
        double s120 = S[1].m[2][0];
        double s121 = S[1].m[2][1];
        double s122 = S[1].m[2][2];

        c1[0] = DerivativeTerm(
            -t0x + t1x,
            (-1 + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
            s000 +
            q0w * q0z * s010 - qperpx * qperpy * s010 +
            qperpw * qperpz * s010 - q0w * q0y * s020 -
            qperpw * qperpy * s020 - qperpx * qperpz * s020 + s100 -
            q0y * q0y * s100 - q0z * q0z * s100 - qperpy * qperpy * s100 -
            qperpz * qperpz * s100 - q0w * q0z * s110 +
            qperpx * qperpy * s110 - qperpw * qperpz * s110 +
            q0w * q0y * s120 + qperpw * qperpy * s120 +
            qperpx * qperpz * s120 +
            q0x * (-(q0y * s010) - q0z * s020 + q0y * s110 + q0z * s120),
            (-1 + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
            s001 +
            q0w * q0z * s011 - qperpx * qperpy * s011 +
            qperpw * qperpz * s011 - q0w * q0y * s021 -
            qperpw * qperpy * s021 - qperpx * qperpz * s021 + s101 -
            q0y * q0y * s101 - q0z * q0z * s101 - qperpy * qperpy * s101 -
            qperpz * qperpz * s101 - q0w * q0z * s111 +
            qperpx * qperpy * s111 - qperpw * qperpz * s111 +
            q0w * q0y * s121 + qperpw * qperpy * s121 +
            qperpx * qperpz * s121 +
            q0x * (-(q0y * s011) - q0z * s021 + q0y * s111 + q0z * s121),
            (-1 + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
            s002 +
            q0w * q0z * s012 - qperpx * qperpy * s012 +
            qperpw * qperpz * s012 - q0w * q0y * s022 -
            qperpw * qperpy * s022 - qperpx * qperpz * s022 + s102 -
            q0y * q0y * s102 - q0z * q0z * s102 - qperpy * qperpy * s102 -
            qperpz * qperpz * s102 - q0w * q0z * s112 +
            qperpx * qperpy * s112 - qperpw * qperpz * s112 +
            q0w * q0y * s122 + qperpw * qperpy * s122 +
            qperpx * qperpz * s122 +
            q0x * (-(q0y * s012) - q0z * s022 + q0y * s112 + q0z * s122));

        c2[0] = DerivativeTerm(
            0.,
            -(qperpy * qperpy * s000) - qperpz * qperpz * s000 +
            qperpx * qperpy * s010 - qperpw * qperpz * s010 +
            qperpw * qperpy * s020 + qperpx * qperpz * s020 +
            q0y * q0y * (s000 - s100) + q0z * q0z * (s000 - s100) +
            qperpy * qperpy * s100 + qperpz * qperpz * s100 -
            qperpx * qperpy * s110 + qperpw * qperpz * s110 -
            qperpw * qperpy * s120 - qperpx * qperpz * s120 +
            2 * q0x * qperpy * s010 * theta -
            2 * q0w * qperpz * s010 * theta +
            2 * q0w * qperpy * s020 * theta +
            2 * q0x * qperpz * s020 * theta +
            q0y *
            (q0x * (-s010 + s110) + q0w * (-s020 + s120) +
                2 * (-2 * qperpy * s000 + qperpx * s010 + qperpw * s020) *
                theta) +
            q0z * (q0w * (s010 - s110) + q0x * (-s020 + s120) -
                2 * (2 * qperpz * s000 + qperpw * s010 - qperpx * s020) *
                theta),
            -(qperpy * qperpy * s001) - qperpz * qperpz * s001 +
            qperpx * qperpy * s011 - qperpw * qperpz * s011 +
            qperpw * qperpy * s021 + qperpx * qperpz * s021 +
            q0y * q0y * (s001 - s101) + q0z * q0z * (s001 - s101) +
            qperpy * qperpy * s101 + qperpz * qperpz * s101 -
            qperpx * qperpy * s111 + qperpw * qperpz * s111 -
            qperpw * qperpy * s121 - qperpx * qperpz * s121 +
            2 * q0x * qperpy * s011 * theta -
            2 * q0w * qperpz * s011 * theta +
            2 * q0w * qperpy * s021 * theta +
            2 * q0x * qperpz * s021 * theta +
            q0y *
            (q0x * (-s011 + s111) + q0w * (-s021 + s121) +
                2 * (-2 * qperpy * s001 + qperpx * s011 + qperpw * s021) *
                theta) +
            q0z * (q0w * (s011 - s111) + q0x * (-s021 + s121) -
                2 * (2 * qperpz * s001 + qperpw * s011 - qperpx * s021) *
                theta),
            -(qperpy * qperpy * s002) - qperpz * qperpz * s002 +
            qperpx * qperpy * s012 - qperpw * qperpz * s012 +
            qperpw * qperpy * s022 + qperpx * qperpz * s022 +
            q0y * q0y * (s002 - s102) + q0z * q0z * (s002 - s102) +
            qperpy * qperpy * s102 + qperpz * qperpz * s102 -
            qperpx * qperpy * s112 + qperpw * qperpz * s112 -
            qperpw * qperpy * s122 - qperpx * qperpz * s122 +
            2 * q0x * qperpy * s012 * theta -
            2 * q0w * qperpz * s012 * theta +
            2 * q0w * qperpy * s022 * theta +
            2 * q0x * qperpz * s022 * theta +
            q0y *
            (q0x * (-s012 + s112) + q0w * (-s022 + s122) +
                2 * (-2 * qperpy * s002 + qperpx * s012 + qperpw * s022) *
                theta) +
            q0z * (q0w * (s012 - s112) + q0x * (-s022 + s122) -
                2 * (2 * qperpz * s002 + qperpw * s012 - qperpx * s022) *
                theta));

        c3[0] = DerivativeTerm(
            0.,
            -2 * (q0x * qperpy * s010 - q0w * qperpz * s010 +
                q0w * qperpy * s020 + q0x * qperpz * s020 -
                q0x * qperpy * s110 + q0w * qperpz * s110 -
                q0w * qperpy * s120 - q0x * qperpz * s120 +
                q0y * (-2 * qperpy * s000 + qperpx * s010 + qperpw * s020 +
                    2 * qperpy * s100 - qperpx * s110 - qperpw * s120) +
                q0z * (-2 * qperpz * s000 - qperpw * s010 + qperpx * s020 +
                    2 * qperpz * s100 + qperpw * s110 - qperpx * s120)) *
            theta,
            -2 * (q0x * qperpy * s011 - q0w * qperpz * s011 +
                q0w * qperpy * s021 + q0x * qperpz * s021 -
                q0x * qperpy * s111 + q0w * qperpz * s111 -
                q0w * qperpy * s121 - q0x * qperpz * s121 +
                q0y * (-2 * qperpy * s001 + qperpx * s011 + qperpw * s021 +
                    2 * qperpy * s101 - qperpx * s111 - qperpw * s121) +
                q0z * (-2 * qperpz * s001 - qperpw * s011 + qperpx * s021 +
                    2 * qperpz * s101 + qperpw * s111 - qperpx * s121)) *
            theta,
            -2 * (q0x * qperpy * s012 - q0w * qperpz * s012 +
                q0w * qperpy * s022 + q0x * qperpz * s022 -
                q0x * qperpy * s112 + q0w * qperpz * s112 -
                q0w * qperpy * s122 - q0x * qperpz * s122 +
                q0y * (-2 * qperpy * s002 + qperpx * s012 + qperpw * s022 +
                    2 * qperpy * s102 - qperpx * s112 - qperpw * s122) +
                q0z * (-2 * qperpz * s002 - qperpw * s012 + qperpx * s022 +
                    2 * qperpz * s102 + qperpw * s112 - qperpx * s122)) *
            theta);

        c4[0] = DerivativeTerm(
            0.,
            -(q0x * qperpy * s010) + q0w * qperpz * s010 - q0w * qperpy * s020 -
            q0x * qperpz * s020 + q0x * qperpy * s110 -
            q0w * qperpz * s110 + q0w * qperpy * s120 +
            q0x * qperpz * s120 + 2 * q0y * q0y * s000 * theta +
            2 * q0z * q0z * s000 * theta -
            2 * qperpy * qperpy * s000 * theta -
            2 * qperpz * qperpz * s000 * theta +
            2 * qperpx * qperpy * s010 * theta -
            2 * qperpw * qperpz * s010 * theta +
            2 * qperpw * qperpy * s020 * theta +
            2 * qperpx * qperpz * s020 * theta +
            q0y * (-(qperpx * s010) - qperpw * s020 +
                2 * qperpy * (s000 - s100) + qperpx * s110 +
                qperpw * s120 - 2 * q0x * s010 * theta -
                2 * q0w * s020 * theta) +
            q0z * (2 * qperpz * s000 + qperpw * s010 - qperpx * s020 -
                2 * qperpz * s100 - qperpw * s110 + qperpx * s120 +
                2 * q0w * s010 * theta - 2 * q0x * s020 * theta),
            -(q0x * qperpy * s011) + q0w * qperpz * s011 - q0w * qperpy * s021 -
            q0x * qperpz * s021 + q0x * qperpy * s111 -
            q0w * qperpz * s111 + q0w * qperpy * s121 +
            q0x * qperpz * s121 + 2 * q0y * q0y * s001 * theta +
            2 * q0z * q0z * s001 * theta -
            2 * qperpy * qperpy * s001 * theta -
            2 * qperpz * qperpz * s001 * theta +
            2 * qperpx * qperpy * s011 * theta -
            2 * qperpw * qperpz * s011 * theta +
            2 * qperpw * qperpy * s021 * theta +
            2 * qperpx * qperpz * s021 * theta +
            q0y * (-(qperpx * s011) - qperpw * s021 +
                2 * qperpy * (s001 - s101) + qperpx * s111 +
                qperpw * s121 - 2 * q0x * s011 * theta -
                2 * q0w * s021 * theta) +
            q0z * (2 * qperpz * s001 + qperpw * s011 - qperpx * s021 -
                2 * qperpz * s101 - qperpw * s111 + qperpx * s121 +
                2 * q0w * s011 * theta - 2 * q0x * s021 * theta),
            -(q0x * qperpy * s012) + q0w * qperpz * s012 - q0w * qperpy * s022 -
            q0x * qperpz * s022 + q0x * qperpy * s112 -
            q0w * qperpz * s112 + q0w * qperpy * s122 +
            q0x * qperpz * s122 + 2 * q0y * q0y * s002 * theta +
            2 * q0z * q0z * s002 * theta -
            2 * qperpy * qperpy * s002 * theta -
            2 * qperpz * qperpz * s002 * theta +
            2 * qperpx * qperpy * s012 * theta -
            2 * qperpw * qperpz * s012 * theta +
            2 * qperpw * qperpy * s022 * theta +
            2 * qperpx * qperpz * s022 * theta +
            q0y * (-(qperpx * s012) - qperpw * s022 +
                2 * qperpy * (s002 - s102) + qperpx * s112 +
                qperpw * s122 - 2 * q0x * s012 * theta -
                2 * q0w * s022 * theta) +
            q0z * (2 * qperpz * s002 + qperpw * s012 - qperpx * s022 -
                2 * qperpz * s102 - qperpw * s112 + qperpx * s122 +
                2 * q0w * s012 * theta - 2 * q0x * s022 * theta));

        c5[0] = DerivativeTerm(
            0.,
            2 * (qperpy * qperpy * s000 + qperpz * qperpz * s000 -
                qperpx * qperpy * s010 + qperpw * qperpz * s010 -
                qperpw * qperpy * s020 - qperpx * qperpz * s020 -
                qperpy * qperpy * s100 - qperpz * qperpz * s100 +
                q0y * q0y * (-s000 + s100) + q0z * q0z * (-s000 + s100) +
                qperpx * qperpy * s110 - qperpw * qperpz * s110 +
                q0y * (q0x * (s010 - s110) + q0w * (s020 - s120)) +
                qperpw * qperpy * s120 + qperpx * qperpz * s120 +
                q0z * (-(q0w * s010) + q0x * s020 + q0w * s110 - q0x * s120)) *
            theta,
            2 * (qperpy * qperpy * s001 + qperpz * qperpz * s001 -
                qperpx * qperpy * s011 + qperpw * qperpz * s011 -
                qperpw * qperpy * s021 - qperpx * qperpz * s021 -
                qperpy * qperpy * s101 - qperpz * qperpz * s101 +
                q0y * q0y * (-s001 + s101) + q0z * q0z * (-s001 + s101) +
                qperpx * qperpy * s111 - qperpw * qperpz * s111 +
                q0y * (q0x * (s011 - s111) + q0w * (s021 - s121)) +
                qperpw * qperpy * s121 + qperpx * qperpz * s121 +
                q0z * (-(q0w * s011) + q0x * s021 + q0w * s111 - q0x * s121)) *
            theta,
            2 * (qperpy * qperpy * s002 + qperpz * qperpz * s002 -
                qperpx * qperpy * s012 + qperpw * qperpz * s012 -
                qperpw * qperpy * s022 - qperpx * qperpz * s022 -
                qperpy * qperpy * s102 - qperpz * qperpz * s102 +
                q0y * q0y * (-s002 + s102) + q0z * q0z * (-s002 + s102) +
                qperpx * qperpy * s112 - qperpw * qperpz * s112 +
                q0y * (q0x * (s012 - s112) + q0w * (s022 - s122)) +
                qperpw * qperpy * s122 + qperpx * qperpz * s122 +
                q0z * (-(q0w * s012) + q0x * s022 + q0w * s112 - q0x * s122)) *
            theta);

        c1[1] = DerivativeTerm(
            -t0y + t1y,
            -(qperpx * qperpy * s000) - qperpw * qperpz * s000 - s010 +
            q0z * q0z * s010 + qperpx * qperpx * s010 +
            qperpz * qperpz * s010 - q0y * q0z * s020 +
            qperpw * qperpx * s020 - qperpy * qperpz * s020 +
            qperpx * qperpy * s100 + qperpw * qperpz * s100 +
            q0w * q0z * (-s000 + s100) + q0x * q0x * (s010 - s110) + s110 -
            q0z * q0z * s110 - qperpx * qperpx * s110 -
            qperpz * qperpz * s110 +
            q0x * (q0y * (-s000 + s100) + q0w * (s020 - s120)) +
            q0y * q0z * s120 - qperpw * qperpx * s120 +
            qperpy * qperpz * s120,
            -(qperpx * qperpy * s001) - qperpw * qperpz * s001 - s011 +
            q0z * q0z * s011 + qperpx * qperpx * s011 +
            qperpz * qperpz * s011 - q0y * q0z * s021 +
            qperpw * qperpx * s021 - qperpy * qperpz * s021 +
            qperpx * qperpy * s101 + qperpw * qperpz * s101 +
            q0w * q0z * (-s001 + s101) + q0x * q0x * (s011 - s111) + s111 -
            q0z * q0z * s111 - qperpx * qperpx * s111 -
            qperpz * qperpz * s111 +
            q0x * (q0y * (-s001 + s101) + q0w * (s021 - s121)) +
            q0y * q0z * s121 - qperpw * qperpx * s121 +
            qperpy * qperpz * s121,
            -(qperpx * qperpy * s002) - qperpw * qperpz * s002 - s012 +
            q0z * q0z * s012 + qperpx * qperpx * s012 +
            qperpz * qperpz * s012 - q0y * q0z * s022 +
            qperpw * qperpx * s022 - qperpy * qperpz * s022 +
            qperpx * qperpy * s102 + qperpw * qperpz * s102 +
            q0w * q0z * (-s002 + s102) + q0x * q0x * (s012 - s112) + s112 -
            q0z * q0z * s112 - qperpx * qperpx * s112 -
            qperpz * qperpz * s112 +
            q0x * (q0y * (-s002 + s102) + q0w * (s022 - s122)) +
            q0y * q0z * s122 - qperpw * qperpx * s122 +
            qperpy * qperpz * s122);

        c2[1] = DerivativeTerm(
            0.,
            qperpx * qperpy * s000 + qperpw * qperpz * s000 + q0z * q0z * s010 -
            qperpx * qperpx * s010 - qperpz * qperpz * s010 -
            q0y * q0z * s020 - qperpw * qperpx * s020 +
            qperpy * qperpz * s020 - qperpx * qperpy * s100 -
            qperpw * qperpz * s100 + q0x * q0x * (s010 - s110) -
            q0z * q0z * s110 + qperpx * qperpx * s110 +
            qperpz * qperpz * s110 + q0y * q0z * s120 +
            qperpw * qperpx * s120 - qperpy * qperpz * s120 +
            2 * q0z * qperpw * s000 * theta +
            2 * q0y * qperpx * s000 * theta -
            4 * q0z * qperpz * s010 * theta +
            2 * q0z * qperpy * s020 * theta +
            2 * q0y * qperpz * s020 * theta +
            q0x * (q0w * s020 + q0y * (-s000 + s100) - q0w * s120 +
                2 * qperpy * s000 * theta - 4 * qperpx * s010 * theta -
                2 * qperpw * s020 * theta) +
            q0w * (-(q0z * s000) + q0z * s100 + 2 * qperpz * s000 * theta -
                2 * qperpx * s020 * theta),
            qperpx * qperpy * s001 + qperpw * qperpz * s001 + q0z * q0z * s011 -
            qperpx * qperpx * s011 - qperpz * qperpz * s011 -
            q0y * q0z * s021 - qperpw * qperpx * s021 +
            qperpy * qperpz * s021 - qperpx * qperpy * s101 -
            qperpw * qperpz * s101 + q0x * q0x * (s011 - s111) -
            q0z * q0z * s111 + qperpx * qperpx * s111 +
            qperpz * qperpz * s111 + q0y * q0z * s121 +
            qperpw * qperpx * s121 - qperpy * qperpz * s121 +
            2 * q0z * qperpw * s001 * theta +
            2 * q0y * qperpx * s001 * theta -
            4 * q0z * qperpz * s011 * theta +
            2 * q0z * qperpy * s021 * theta +
            2 * q0y * qperpz * s021 * theta +
            q0x * (q0w * s021 + q0y * (-s001 + s101) - q0w * s121 +
                2 * qperpy * s001 * theta - 4 * qperpx * s011 * theta -
                2 * qperpw * s021 * theta) +
            q0w * (-(q0z * s001) + q0z * s101 + 2 * qperpz * s001 * theta -
                2 * qperpx * s021 * theta),
            qperpx * qperpy * s002 + qperpw * qperpz * s002 + q0z * q0z * s012 -
            qperpx * qperpx * s012 - qperpz * qperpz * s012 -
            q0y * q0z * s022 - qperpw * qperpx * s022 +
            qperpy * qperpz * s022 - qperpx * qperpy * s102 -
            qperpw * qperpz * s102 + q0x * q0x * (s012 - s112) -
            q0z * q0z * s112 + qperpx * qperpx * s112 +
            qperpz * qperpz * s112 + q0y * q0z * s122 +
            qperpw * qperpx * s122 - qperpy * qperpz * s122 +
            2 * q0z * qperpw * s002 * theta +
            2 * q0y * qperpx * s002 * theta -
            4 * q0z * qperpz * s012 * theta +
            2 * q0z * qperpy * s022 * theta +
            2 * q0y * qperpz * s022 * theta +
            q0x * (q0w * s022 + q0y * (-s002 + s102) - q0w * s122 +
                2 * qperpy * s002 * theta - 4 * qperpx * s012 * theta -
                2 * qperpw * s022 * theta) +
            q0w * (-(q0z * s002) + q0z * s102 + 2 * qperpz * s002 * theta -
                2 * qperpx * s022 * theta));

        c3[1] = DerivativeTerm(
            0., 2 * (-(q0x * qperpy * s000) - q0w * qperpz * s000 +
                2 * q0x * qperpx * s010 + q0x * qperpw * s020 +
                q0w * qperpx * s020 + q0x * qperpy * s100 +
                q0w * qperpz * s100 - 2 * q0x * qperpx * s110 -
                q0x * qperpw * s120 - q0w * qperpx * s120 +
                q0z * (2 * qperpz * s010 - qperpy * s020 +
                    qperpw * (-s000 + s100) - 2 * qperpz * s110 +
                    qperpy * s120) +
                q0y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 +
                    qperpz * s120)) *
            theta,
            2 * (-(q0x * qperpy * s001) - q0w * qperpz * s001 +
                2 * q0x * qperpx * s011 + q0x * qperpw * s021 +
                q0w * qperpx * s021 + q0x * qperpy * s101 +
                q0w * qperpz * s101 - 2 * q0x * qperpx * s111 -
                q0x * qperpw * s121 - q0w * qperpx * s121 +
                q0z * (2 * qperpz * s011 - qperpy * s021 +
                    qperpw * (-s001 + s101) - 2 * qperpz * s111 +
                    qperpy * s121) +
                q0y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 +
                    qperpz * s121)) *
            theta,
            2 * (-(q0x * qperpy * s002) - q0w * qperpz * s002 +
                2 * q0x * qperpx * s012 + q0x * qperpw * s022 +
                q0w * qperpx * s022 + q0x * qperpy * s102 +
                q0w * qperpz * s102 - 2 * q0x * qperpx * s112 -
                q0x * qperpw * s122 - q0w * qperpx * s122 +
                q0z * (2 * qperpz * s012 - qperpy * s022 +
                    qperpw * (-s002 + s102) - 2 * qperpz * s112 +
                    qperpy * s122) +
                q0y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 +
                    qperpz * s122)) *
            theta);

        c4[1] = DerivativeTerm(
            0.,
            -(q0x * qperpy * s000) - q0w * qperpz * s000 +
            2 * q0x * qperpx * s010 + q0x * qperpw * s020 +
            q0w * qperpx * s020 + q0x * qperpy * s100 +
            q0w * qperpz * s100 - 2 * q0x * qperpx * s110 -
            q0x * qperpw * s120 - q0w * qperpx * s120 +
            2 * qperpx * qperpy * s000 * theta +
            2 * qperpw * qperpz * s000 * theta +
            2 * q0x * q0x * s010 * theta + 2 * q0z * q0z * s010 * theta -
            2 * qperpx * qperpx * s010 * theta -
            2 * qperpz * qperpz * s010 * theta +
            2 * q0w * q0x * s020 * theta -
            2 * qperpw * qperpx * s020 * theta +
            2 * qperpy * qperpz * s020 * theta +
            q0y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 +
                qperpz * s120 - 2 * q0x * s000 * theta) +
            q0z * (2 * qperpz * s010 - qperpy * s020 +
                qperpw * (-s000 + s100) - 2 * qperpz * s110 +
                qperpy * s120 - 2 * q0w * s000 * theta -
                2 * q0y * s020 * theta),
            -(q0x * qperpy * s001) - q0w * qperpz * s001 +
            2 * q0x * qperpx * s011 + q0x * qperpw * s021 +
            q0w * qperpx * s021 + q0x * qperpy * s101 +
            q0w * qperpz * s101 - 2 * q0x * qperpx * s111 -
            q0x * qperpw * s121 - q0w * qperpx * s121 +
            2 * qperpx * qperpy * s001 * theta +
            2 * qperpw * qperpz * s001 * theta +
            2 * q0x * q0x * s011 * theta + 2 * q0z * q0z * s011 * theta -
            2 * qperpx * qperpx * s011 * theta -
            2 * qperpz * qperpz * s011 * theta +
            2 * q0w * q0x * s021 * theta -
            2 * qperpw * qperpx * s021 * theta +
            2 * qperpy * qperpz * s021 * theta +
            q0y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 +
                qperpz * s121 - 2 * q0x * s001 * theta) +
            q0z * (2 * qperpz * s011 - qperpy * s021 +
                qperpw * (-s001 + s101) - 2 * qperpz * s111 +
                qperpy * s121 - 2 * q0w * s001 * theta -
                2 * q0y * s021 * theta),
            -(q0x * qperpy * s002) - q0w * qperpz * s002 +
            2 * q0x * qperpx * s012 + q0x * qperpw * s022 +
            q0w * qperpx * s022 + q0x * qperpy * s102 +
            q0w * qperpz * s102 - 2 * q0x * qperpx * s112 -
            q0x * qperpw * s122 - q0w * qperpx * s122 +
            2 * qperpx * qperpy * s002 * theta +
            2 * qperpw * qperpz * s002 * theta +
            2 * q0x * q0x * s012 * theta + 2 * q0z * q0z * s012 * theta -
            2 * qperpx * qperpx * s012 * theta -
            2 * qperpz * qperpz * s012 * theta +
            2 * q0w * q0x * s022 * theta -
            2 * qperpw * qperpx * s022 * theta +
            2 * qperpy * qperpz * s022 * theta +
            q0y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 +
                qperpz * s122 - 2 * q0x * s002 * theta) +
            q0z * (2 * qperpz * s012 - qperpy * s022 +
                qperpw * (-s002 + s102) - 2 * qperpz * s112 +
                qperpy * s122 - 2 * q0w * s002 * theta -
                2 * q0y * s022 * theta));

        c5[1] = DerivativeTerm(
            0., -2 * (qperpx * qperpy * s000 + qperpw * qperpz * s000 +
                q0z * q0z * s010 - qperpx * qperpx * s010 -
                qperpz * qperpz * s010 - q0y * q0z * s020 -
                qperpw * qperpx * s020 + qperpy * qperpz * s020 -
                qperpx * qperpy * s100 - qperpw * qperpz * s100 +
                q0w * q0z * (-s000 + s100) + q0x * q0x * (s010 - s110) -
                q0z * q0z * s110 + qperpx * qperpx * s110 +
                qperpz * qperpz * s110 +
                q0x * (q0y * (-s000 + s100) + q0w * (s020 - s120)) +
                q0y * q0z * s120 + qperpw * qperpx * s120 -
                qperpy * qperpz * s120) *
            theta,
            -2 * (qperpx * qperpy * s001 + qperpw * qperpz * s001 +
                q0z * q0z * s011 - qperpx * qperpx * s011 -
                qperpz * qperpz * s011 - q0y * q0z * s021 -
                qperpw * qperpx * s021 + qperpy * qperpz * s021 -
                qperpx * qperpy * s101 - qperpw * qperpz * s101 +
                q0w * q0z * (-s001 + s101) + q0x * q0x * (s011 - s111) -
                q0z * q0z * s111 + qperpx * qperpx * s111 +
                qperpz * qperpz * s111 +
                q0x * (q0y * (-s001 + s101) + q0w * (s021 - s121)) +
                q0y * q0z * s121 + qperpw * qperpx * s121 -
                qperpy * qperpz * s121) *
            theta,
            -2 * (qperpx * qperpy * s002 + qperpw * qperpz * s002 +
                q0z * q0z * s012 - qperpx * qperpx * s012 -
                qperpz * qperpz * s012 - q0y * q0z * s022 -
                qperpw * qperpx * s022 + qperpy * qperpz * s022 -
                qperpx * qperpy * s102 - qperpw * qperpz * s102 +
                q0w * q0z * (-s002 + s102) + q0x * q0x * (s012 - s112) -
                q0z * q0z * s112 + qperpx * qperpx * s112 +
                qperpz * qperpz * s112 +
                q0x * (q0y * (-s002 + s102) + q0w * (s022 - s122)) +
                q0y * q0z * s122 + qperpw * qperpx * s122 -
                qperpy * qperpz * s122) *
            theta);

        c1[2] = DerivativeTerm(
            -t0z + t1z, (qperpw * qperpy * s000 - qperpx * qperpz * s000 -
                q0y * q0z * s010 - qperpw * qperpx * s010 -
                qperpy * qperpz * s010 - s020 + q0y * q0y * s020 +
                qperpx * qperpx * s020 + qperpy * qperpy * s020 -
                qperpw * qperpy * s100 + qperpx * qperpz * s100 +
                q0x * q0z * (-s000 + s100) + q0y * q0z * s110 +
                qperpw * qperpx * s110 + qperpy * qperpz * s110 +
                q0w * (q0y * (s000 - s100) + q0x * (-s010 + s110)) +
                q0x * q0x * (s020 - s120) + s120 - q0y * q0y * s120 -
                qperpx * qperpx * s120 - qperpy * qperpy * s120),
            (qperpw * qperpy * s001 - qperpx * qperpz * s001 -
                q0y * q0z * s011 - qperpw * qperpx * s011 -
                qperpy * qperpz * s011 - s021 + q0y * q0y * s021 +
                qperpx * qperpx * s021 + qperpy * qperpy * s021 -
                qperpw * qperpy * s101 + qperpx * qperpz * s101 +
                q0x * q0z * (-s001 + s101) + q0y * q0z * s111 +
                qperpw * qperpx * s111 + qperpy * qperpz * s111 +
                q0w * (q0y * (s001 - s101) + q0x * (-s011 + s111)) +
                q0x * q0x * (s021 - s121) + s121 - q0y * q0y * s121 -
                qperpx * qperpx * s121 - qperpy * qperpy * s121),
            (qperpw * qperpy * s002 - qperpx * qperpz * s002 -
                q0y * q0z * s012 - qperpw * qperpx * s012 -
                qperpy * qperpz * s012 - s022 + q0y * q0y * s022 +
                qperpx * qperpx * s022 + qperpy * qperpy * s022 -
                qperpw * qperpy * s102 + qperpx * qperpz * s102 +
                q0x * q0z * (-s002 + s102) + q0y * q0z * s112 +
                qperpw * qperpx * s112 + qperpy * qperpz * s112 +
                q0w * (q0y * (s002 - s102) + q0x * (-s012 + s112)) +
                q0x * q0x * (s022 - s122) + s122 - q0y * q0y * s122 -
                qperpx * qperpx * s122 - qperpy * qperpy * s122));

        c2[2] = DerivativeTerm(
            0.,
            (q0w * q0y * s000 - q0x * q0z * s000 - qperpw * qperpy * s000 +
                qperpx * qperpz * s000 - q0w * q0x * s010 - q0y * q0z * s010 +
                qperpw * qperpx * s010 + qperpy * qperpz * s010 +
                q0x * q0x * s020 + q0y * q0y * s020 - qperpx * qperpx * s020 -
                qperpy * qperpy * s020 - q0w * q0y * s100 + q0x * q0z * s100 +
                qperpw * qperpy * s100 - qperpx * qperpz * s100 +
                q0w * q0x * s110 + q0y * q0z * s110 - qperpw * qperpx * s110 -
                qperpy * qperpz * s110 - q0x * q0x * s120 - q0y * q0y * s120 +
                qperpx * qperpx * s120 + qperpy * qperpy * s120 -
                2 * q0y * qperpw * s000 * theta + 2 * q0z * qperpx * s000 * theta -
                2 * q0w * qperpy * s000 * theta + 2 * q0x * qperpz * s000 * theta +
                2 * q0x * qperpw * s010 * theta + 2 * q0w * qperpx * s010 * theta +
                2 * q0z * qperpy * s010 * theta + 2 * q0y * qperpz * s010 * theta -
                4 * q0x * qperpx * s020 * theta - 4 * q0y * qperpy * s020 * theta),
            (q0w * q0y * s001 - q0x * q0z * s001 - qperpw * qperpy * s001 +
                qperpx * qperpz * s001 - q0w * q0x * s011 - q0y * q0z * s011 +
                qperpw * qperpx * s011 + qperpy * qperpz * s011 +
                q0x * q0x * s021 + q0y * q0y * s021 - qperpx * qperpx * s021 -
                qperpy * qperpy * s021 - q0w * q0y * s101 + q0x * q0z * s101 +
                qperpw * qperpy * s101 - qperpx * qperpz * s101 +
                q0w * q0x * s111 + q0y * q0z * s111 - qperpw * qperpx * s111 -
                qperpy * qperpz * s111 - q0x * q0x * s121 - q0y * q0y * s121 +
                qperpx * qperpx * s121 + qperpy * qperpy * s121 -
                2 * q0y * qperpw * s001 * theta + 2 * q0z * qperpx * s001 * theta -
                2 * q0w * qperpy * s001 * theta + 2 * q0x * qperpz * s001 * theta +
                2 * q0x * qperpw * s011 * theta + 2 * q0w * qperpx * s011 * theta +
                2 * q0z * qperpy * s011 * theta + 2 * q0y * qperpz * s011 * theta -
                4 * q0x * qperpx * s021 * theta - 4 * q0y * qperpy * s021 * theta),
            (q0w * q0y * s002 - q0x * q0z * s002 - qperpw * qperpy * s002 +
                qperpx * qperpz * s002 - q0w * q0x * s012 - q0y * q0z * s012 +
                qperpw * qperpx * s012 + qperpy * qperpz * s012 +
                q0x * q0x * s022 + q0y * q0y * s022 - qperpx * qperpx * s022 -
                qperpy * qperpy * s022 - q0w * q0y * s102 + q0x * q0z * s102 +
                qperpw * qperpy * s102 - qperpx * qperpz * s102 +
                q0w * q0x * s112 + q0y * q0z * s112 - qperpw * qperpx * s112 -
                qperpy * qperpz * s112 - q0x * q0x * s122 - q0y * q0y * s122 +
                qperpx * qperpx * s122 + qperpy * qperpy * s122 -
                2 * q0y * qperpw * s002 * theta + 2 * q0z * qperpx * s002 * theta -
                2 * q0w * qperpy * s002 * theta + 2 * q0x * qperpz * s002 * theta +
                2 * q0x * qperpw * s012 * theta + 2 * q0w * qperpx * s012 * theta +
                2 * q0z * qperpy * s012 * theta + 2 * q0y * qperpz * s012 * theta -
                4 * q0x * qperpx * s022 * theta -
                4 * q0y * qperpy * s022 * theta));

        c3[2] = DerivativeTerm(
            0., -2 * (-(q0w * qperpy * s000) + q0x * qperpz * s000 +
                q0x * qperpw * s010 + q0w * qperpx * s010 -
                2 * q0x * qperpx * s020 + q0w * qperpy * s100 -
                q0x * qperpz * s100 - q0x * qperpw * s110 -
                q0w * qperpx * s110 +
                q0z * (qperpx * s000 + qperpy * s010 - qperpx * s100 -
                    qperpy * s110) +
                2 * q0x * qperpx * s120 +
                q0y * (qperpz * s010 - 2 * qperpy * s020 +
                    qperpw * (-s000 + s100) - qperpz * s110 +
                    2 * qperpy * s120)) *
            theta,
            -2 * (-(q0w * qperpy * s001) + q0x * qperpz * s001 +
                q0x * qperpw * s011 + q0w * qperpx * s011 -
                2 * q0x * qperpx * s021 + q0w * qperpy * s101 -
                q0x * qperpz * s101 - q0x * qperpw * s111 -
                q0w * qperpx * s111 +
                q0z * (qperpx * s001 + qperpy * s011 - qperpx * s101 -
                    qperpy * s111) +
                2 * q0x * qperpx * s121 +
                q0y * (qperpz * s011 - 2 * qperpy * s021 +
                    qperpw * (-s001 + s101) - qperpz * s111 +
                    2 * qperpy * s121)) *
            theta,
            -2 * (-(q0w * qperpy * s002) + q0x * qperpz * s002 +
                q0x * qperpw * s012 + q0w * qperpx * s012 -
                2 * q0x * qperpx * s022 + q0w * qperpy * s102 -
                q0x * qperpz * s102 - q0x * qperpw * s112 -
                q0w * qperpx * s112 +
                q0z * (qperpx * s002 + qperpy * s012 - qperpx * s102 -
                    qperpy * s112) +
                2 * q0x * qperpx * s122 +
                q0y * (qperpz * s012 - 2 * qperpy * s022 +
                    qperpw * (-s002 + s102) - qperpz * s112 +
                    2 * qperpy * s122)) *
            theta);

        c4[2] = DerivativeTerm(
            0.,
            q0w * qperpy * s000 - q0x * qperpz * s000 - q0x * qperpw * s010 -
            q0w * qperpx * s010 + 2 * q0x * qperpx * s020 -
            q0w * qperpy * s100 + q0x * qperpz * s100 +
            q0x * qperpw * s110 + q0w * qperpx * s110 -
            2 * q0x * qperpx * s120 - 2 * qperpw * qperpy * s000 * theta +
            2 * qperpx * qperpz * s000 * theta -
            2 * q0w * q0x * s010 * theta +
            2 * qperpw * qperpx * s010 * theta +
            2 * qperpy * qperpz * s010 * theta +
            2 * q0x * q0x * s020 * theta + 2 * q0y * q0y * s020 * theta -
            2 * qperpx * qperpx * s020 * theta -
            2 * qperpy * qperpy * s020 * theta +
            q0z * (-(qperpx * s000) - qperpy * s010 + qperpx * s100 +
                qperpy * s110 - 2 * q0x * s000 * theta) +
            q0y * (-(qperpz * s010) + 2 * qperpy * s020 +
                qperpw * (s000 - s100) + qperpz * s110 -
                2 * qperpy * s120 + 2 * q0w * s000 * theta -
                2 * q0z * s010 * theta),
            q0w * qperpy * s001 - q0x * qperpz * s001 - q0x * qperpw * s011 -
            q0w * qperpx * s011 + 2 * q0x * qperpx * s021 -
            q0w * qperpy * s101 + q0x * qperpz * s101 +
            q0x * qperpw * s111 + q0w * qperpx * s111 -
            2 * q0x * qperpx * s121 - 2 * qperpw * qperpy * s001 * theta +
            2 * qperpx * qperpz * s001 * theta -
            2 * q0w * q0x * s011 * theta +
            2 * qperpw * qperpx * s011 * theta +
            2 * qperpy * qperpz * s011 * theta +
            2 * q0x * q0x * s021 * theta + 2 * q0y * q0y * s021 * theta -
            2 * qperpx * qperpx * s021 * theta -
            2 * qperpy * qperpy * s021 * theta +
            q0z * (-(qperpx * s001) - qperpy * s011 + qperpx * s101 +
                qperpy * s111 - 2 * q0x * s001 * theta) +
            q0y * (-(qperpz * s011) + 2 * qperpy * s021 +
                qperpw * (s001 - s101) + qperpz * s111 -
                2 * qperpy * s121 + 2 * q0w * s001 * theta -
                2 * q0z * s011 * theta),
            q0w * qperpy * s002 - q0x * qperpz * s002 - q0x * qperpw * s012 -
            q0w * qperpx * s012 + 2 * q0x * qperpx * s022 -
            q0w * qperpy * s102 + q0x * qperpz * s102 +
            q0x * qperpw * s112 + q0w * qperpx * s112 -
            2 * q0x * qperpx * s122 - 2 * qperpw * qperpy * s002 * theta +
            2 * qperpx * qperpz * s002 * theta -
            2 * q0w * q0x * s012 * theta +
            2 * qperpw * qperpx * s012 * theta +
            2 * qperpy * qperpz * s012 * theta +
            2 * q0x * q0x * s022 * theta + 2 * q0y * q0y * s022 * theta -
            2 * qperpx * qperpx * s022 * theta -
            2 * qperpy * qperpy * s022 * theta +
            q0z * (-(qperpx * s002) - qperpy * s012 + qperpx * s102 +
                qperpy * s112 - 2 * q0x * s002 * theta) +
            q0y * (-(qperpz * s012) + 2 * qperpy * s022 +
                qperpw * (s002 - s102) + qperpz * s112 -
                2 * qperpy * s122 + 2 * q0w * s002 * theta -
                2 * q0z * s012 * theta));

        c5[2] = DerivativeTerm(
            0., 2 * (qperpw * qperpy * s000 - qperpx * qperpz * s000 +
                q0y * q0z * s010 - qperpw * qperpx * s010 -
                qperpy * qperpz * s010 - q0y * q0y * s020 +
                qperpx * qperpx * s020 + qperpy * qperpy * s020 +
                q0x * q0z * (s000 - s100) - qperpw * qperpy * s100 +
                qperpx * qperpz * s100 +
                q0w * (q0y * (-s000 + s100) + q0x * (s010 - s110)) -
                q0y * q0z * s110 + qperpw * qperpx * s110 +
                qperpy * qperpz * s110 + q0y * q0y * s120 -
                qperpx * qperpx * s120 - qperpy * qperpy * s120 +
                q0x * q0x * (-s020 + s120)) *
            theta,
            2 * (qperpw * qperpy * s001 - qperpx * qperpz * s001 +
                q0y * q0z * s011 - qperpw * qperpx * s011 -
                qperpy * qperpz * s011 - q0y * q0y * s021 +
                qperpx * qperpx * s021 + qperpy * qperpy * s021 +
                q0x * q0z * (s001 - s101) - qperpw * qperpy * s101 +
                qperpx * qperpz * s101 +
                q0w * (q0y * (-s001 + s101) + q0x * (s011 - s111)) -
                q0y * q0z * s111 + qperpw * qperpx * s111 +
                qperpy * qperpz * s111 + q0y * q0y * s121 -
                qperpx * qperpx * s121 - qperpy * qperpy * s121 +
                q0x * q0x * (-s021 + s121)) *
            theta,
            2 * (qperpw * qperpy * s002 - qperpx * qperpz * s002 +
                q0y * q0z * s012 - qperpw * qperpx * s012 -
                qperpy * qperpz * s012 - q0y * q0y * s022 +
                qperpx * qperpx * s022 + qperpy * qperpy * s022 +
                q0x * q0z * (s002 - s102) - qperpw * qperpy * s102 +
                qperpx * qperpz * s102 +
                q0w * (q0y * (-s002 + s102) + q0x * (s012 - s112)) -
                q0y * q0z * s112 + qperpw * qperpx * s112 +
                qperpy * qperpz * s112 + q0y * q0y * s122 -
                qperpx * qperpx * s122 - qperpy * qperpy * s122 +
                q0x * q0x * (-s022 + s122)) *
            theta);
    }
}

void AnimatedTransform::Decompose(const Matrix4x4& m, Vector3f* T,
    Quaternion* Rquat, Matrix4x4* S) {
    // Extract translation _T_ from transformation matrix
    T->x = m.m[0][3];
    T->y = m.m[1][3];
    T->z = m.m[2][3];

    // Compute new transformation matrix _M_ without translation
    Matrix4x4 M = m;
    for (int i = 0; i < 3; ++i) M.m[i][3] = M.m[3][i] = 0.f;
    M.m[3][3] = 1.f;

    // Extract rotation _R_ from transformation matrix
    double norm;
    int count = 0;
    Matrix4x4 R = M;
    do {
        // Compute next matrix _Rnext_ in series
        Matrix4x4 Rnext;
        Matrix4x4 Rit = Inverse(Transpose(R));
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                Rnext.m[i][j] = 0.5f * (R.m[i][j] + Rit.m[i][j]);

        // Compute norm of difference between _R_ and _Rnext_
        norm = 0;
        for (int i = 0; i < 3; ++i) {
            double n = std::abs(R.m[i][0] - Rnext.m[i][0]) +
                std::abs(R.m[i][1] - Rnext.m[i][1]) +
                std::abs(R.m[i][2] - Rnext.m[i][2]);
            norm = std::max(norm, n);
        }
        R = Rnext;
    } while (++count < 100 && norm > .0001);
    // XXX TODO FIXME deal with flip...
    *Rquat = Quaternion(R);

    // Compute scale _S_ using rotation and original matrix
    *S = Matrix4x4::Mul(Inverse(R), M);
}

void AnimatedTransform::Interpolate(double time, Transform* t) const {
    // Handle boundary conditions for matrix interpolation
    if (!actuallyAnimated || time <= startTime) {
        *t = *startTransform;
        return;
    }
    if (time >= endTime) {
        *t = *endTransform;
        return;
    }
    double dt = (time - startTime) / (endTime - startTime);
    // Interpolate translation at _dt_
    Vector3f trans = (1 - dt) * T[0] + dt * T[1];

    // Interpolate rotation at _dt_
    Quaternion rotate = Slerp(dt, R[0], R[1]);

    // Interpolate scale at _dt_
    Matrix4x4 scale;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            scale.m[i][j] = check::Lerp(dt, S[0].m[i][j], S[1].m[i][j]);

    // Compute interpolated matrix as product of interpolated components
    *t = Translate(trans) * rotate.ToTransform() * Transform(scale);
}

Ray AnimatedTransform::operator()(const Ray& r) const {
    if (!actuallyAnimated || r.time <= startTime)
        return (*startTransform)(r);
    else if (r.time >= endTime)
        return (*endTransform)(r);
    else {
        Transform t;
        Interpolate(r.time, &t);
        return t(r);
    }
}

RayDifferential AnimatedTransform::operator()(const RayDifferential& r) const {
    if (!actuallyAnimated || r.time <= startTime)
        return (*startTransform)(r);
    else if (r.time >= endTime)
        return (*endTransform)(r);
    else {
        Transform t;
        Interpolate(r.time, &t);
        return t(r);
    }
}

Point3f AnimatedTransform::operator()(double time, const Point3f& p) const {
    if (!actuallyAnimated || time <= startTime)
        return (*startTransform)(p);
    else if (time >= endTime)
        return (*endTransform)(p);
    Transform t;
    Interpolate(time, &t);
    return t(p);
}

Vector3f AnimatedTransform::operator()(double time, const Vector3f& v) const {
    if (!actuallyAnimated || time <= startTime)
        return (*startTransform)(v);
    else if (time >= endTime)
        return (*endTransform)(v);
    Transform t;
    Interpolate(time, &t);
    return t(v);
}

Bounds3f AnimatedTransform::MotionBounds(const Bounds3f& b) const {
    if (!actuallyAnimated) return (*startTransform)(b);
    if (hasRotation == false)
        return Union((*startTransform)(b), (*endTransform)(b));
    // Return motion bounds accounting for animated rotation
    Bounds3f bounds;
    for (int corner = 0; corner < 8; ++corner)
        bounds = Union(bounds, BoundPointMotion(b.Corner(corner)));
    return bounds;
}*/

inline Interval Sin(const Interval& i) {
    check::CHECK_GE(i.low, 0);
    check::CHECK_LE(i.high, 2.0001 * PI);
    double sinLow = std::sin(i.low), sinHigh = std::sin(i.high);
    if (sinLow > sinHigh) std::swap(sinLow, sinHigh);
    if (i.low < PI / 2 && i.high > PI / 2) sinHigh = 1.;
    if (i.low < (3.f / 2.f) * PI && i.high >(3.f / 2.f) * PI) sinLow = -1.;
    return Interval(sinLow, sinHigh);
}

inline Interval Cos(const Interval& i) {
    check::CHECK_GE(i.low, 0);
    check::CHECK_LE(i.high, 2.0001 * PI);
    double cosLow = std::cos(i.low), cosHigh = std::cos(i.high);
    if (cosLow > cosHigh) std::swap(cosLow, cosHigh);
    if (i.low < PI && i.high > PI) cosLow = -1.;
    return Interval(cosLow, cosHigh);
}


void IntervalFindZeros(double c1, double c2, double c3, double c4, double c5,
    double theta, Interval tInterval, double* zeros,
    int* zeroCount, int depth = 8) {
    // Evaluate motion derivative in interval form, return if no zeros
    Interval range = Interval(c1) +
        (Interval(c2) + Interval(c3) * tInterval) *
        Cos(Interval(2 * theta) * tInterval) +
        (Interval(c4) + Interval(c5) * tInterval) *
        Sin(Interval(2 * theta) * tInterval);
    if (range.low > 0. || range.high < 0. || range.low == range.high) return;
    if (depth > 0) {
        // Split _tInterval_ and check both resulting intervals
        double mid = (tInterval.low + tInterval.high) * 0.5f;
        IntervalFindZeros(c1, c2, c3, c4, c5, theta,
            Interval(tInterval.low, mid), zeros, zeroCount,
            depth - 1);
        IntervalFindZeros(c1, c2, c3, c4, c5, theta,
            Interval(mid, tInterval.high), zeros, zeroCount,
            depth - 1);
    }
    else {
        // Use Newton's method to refine zero
        double tNewton = (tInterval.low + tInterval.high) * 0.5f;
        for (int i = 0; i < 4; ++i) {
            double fNewton =
                c1 + (c2 + c3 * tNewton) * std::cos(2.f * theta * tNewton) +
                (c4 + c5 * tNewton) * std::sin(2.f * theta * tNewton);
            double fPrimeNewton = (c3 + 2 * (c4 + c5 * tNewton) * theta) *
                std::cos(2.f * tNewton * theta) +
                (c5 - 2 * (c2 + c3 * tNewton) * theta) *
                std::sin(2.f * tNewton * theta);
            if (fNewton == 0 || fPrimeNewton == 0) break;
            tNewton = tNewton - fNewton / fPrimeNewton;
        }
        if (tNewton >= tInterval.low - 1e-3f &&
            tNewton < tInterval.high + 1e-3f) {
            zeros[*zeroCount] = tNewton;
            (*zeroCount)++;
        }
    }
}


/*Bounds3f AnimatedTransform::BoundPointMotion(const Point3f& p) const {
    if (!actuallyAnimated) return Bounds3f((*startTransform)(p));
    Bounds3f bounds((*startTransform)(p), (*endTransform)(p));
    double cosTheta = Dot(R[0], R[1]);
    double theta = std::acos(check::Clamp(cosTheta, -1, 1));
    for (int c = 0; c < 3; ++c) {
        // Find any motion derivative zeros for the component _c_
        double zeros[8];
        int nZeros = 0;
        IntervalFindZeros(c1[c].Eval(p), c2[c].Eval(p), c3[c].Eval(p),
            c4[c].Eval(p), c5[c].Eval(p), theta, Interval(0., 1.),
            zeros, &nZeros);
        check::CHECK_LE(nZeros, sizeof(zeros) / sizeof(zeros[0]));

        // Expand bounding box for any motion derivative zeros found
        for (int i = 0; i < nZeros; ++i) {
            Point3f pz = (*this)(check::Lerp(zeros[i], startTime, endTime), p);
            bounds = Union(bounds, pz);
        }
    }
    return bounds;
}*/


