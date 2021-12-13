#ifndef SPHERE_H
#define SPHERE_H

//#include "hittable.h"
#include "onb.h"
#include "efloat.h"
#include "primitive.h"
#include "material.h"
#include "transform.h"

/*class sphere : public Primitive {
public:
    sphere() {}
    sphere(Point3f cen, double r, shared_ptr<material> m)
        : center(cen), radius(r), mat_ptr(m) {};
    sphere(Point3f cen, double r, shared_ptr<material> m,double phi,double thetaM,double thetam)
        : center(cen), radius(r), mat_ptr(m),phiMax(phi),thetaMax(thetaM),thetaMin(thetam) {};

    virtual bool Intersect(
        const ray& r, double t_min, double t_max, SurfaceInteraction& inter_record) const override;
    virtual bool IntersectP(const ray& r, double t_min, double t_max) const override;
    virtual aabb bounding_box() const override;
    double pdf_value(const Point3f& o, const Vector3f& v) const;
    Vector3f random(const Point3f& o) const;

public:
    Point3f center;
    double radius;
    shared_ptr<material> mat_ptr;
    double phiMax=2*PI ;
    double thetaMax=PI , thetaMin=0;

private:
    static void get_sphere_uv(const Point3f& p, double& u, double& v); 
};*/

class sphere : public Primitive {
public:
    /*sphere() {}
    sphere(Point3f cen, double r, shared_ptr<material> m)
        : center(cen), radius(r), mat_ptr(m) {};
    sphere(Point3f cen, double r, shared_ptr<material> m, double phi, double thetaM, double thetam)
        : center(cen), radius(r), mat_ptr(m), phiMax(phi), thetaMax(thetaM), thetaMin(thetam) {};*/
    sphere(Point3f cen, double r, shared_ptr<material> m)
        : center(cen), radius(r), mat_ptr(m),zMax(r),zMin(-r), reverseOrientation(false){};
    sphere(const Transform &OToWorld, const Transform &WToObject,
        bool reverse, Float radius, Float zMin, Float zMax,
        Float phiMax, shared_ptr<material> m)
        : ObjectToWorld(OToWorld),
        WorldToObject(WToObject),
        reverseOrientation(reverse),
        radius(radius),
        zMin(check::Clamp(std::min(zMin, zMax), -radius, radius)),
        zMax(check::Clamp(std::max(zMin, zMax), -radius, radius)),
        thetaMin(std::acos(check::Clamp(std::min(zMin, zMax) / radius, -1, 1))),
        thetaMax(std::acos(check::Clamp(std::max(zMin, zMax) / radius, -1, 1))),
        phiMax(Radians(check::Clamp(phiMax, 0, 360))),
        mat_ptr(m),
        center(ConvertPTrans(OToWorld, Point3f(0, 0, 0)) ){};

    virtual bool Intersect(
        const ray& r, double t_min, double t_max, SurfaceInteraction& inter_record) const override;
    virtual bool IntersectP(const ray& r, double t_min, double t_max) const override;
    virtual aabb bounding_box() const override;
    double pdf_value(const Point3f& o, const Vector3f& v) const;
    Vector3f random(const Point3f& o) const;

public:
    Point3f center;
    double radius;
    shared_ptr<material> mat_ptr;
    double phiMax = 2 * PI;
    double thetaMax = PI, thetaMin = 0;

    //const Float radius;
    const Float zMin, zMax;
    const Transform ObjectToWorld,  WorldToObject;
    const bool reverseOrientation;
    //const bool transformSwapsHandedness;
    // const Float thetaMin, thetaMax, phiMax;

private:
    static void get_sphere_uv(const Point3f& p, double& u, double& v);
};

inline Vector3f random_to_sphere(double radius, double distance_squared) {
    auto r1 = random_double();
    auto r2 = random_double();
    auto z = 1 + r2 * (sqrt(1 - radius * radius / distance_squared) - 1);

    auto phi = 2 * PI * r1;
    auto x = cos(phi) * sqrt(1 - z * z);
    auto y = sin(phi) * sqrt(1 - z * z);

    return Vector3f(x, y, z);
}

#endif
