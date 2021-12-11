#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
#include "shape.h"

#include "efloat.h"

class sphere : public hittable {
public:
    sphere() {}
    sphere(point3 cen, Float r, shared_ptr<material> m)
        : center(cen), radius(r), mat_ptr(m) {};

    virtual bool hit(
        const ray& r, Float t_min, Float t_max, hit_record& rec) const override;
    virtual bool bounding_box(Float time0, Float time1, aabb& output_box) const override;
public:
    point3 center;
    Float radius;
    shared_ptr<material> mat_ptr;
private:
    static void get_sphere_uv(const point3& p, Float& u, Float& v) {
        // p: a given point on the sphere of radius one, centered at the origin.
        // u: returned value [0,1] of angle around the Y axis from X=-1.
        // v: returned value [0,1] of angle from Y=-1 to Y=+1.
        //     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
        //     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
        //     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>

        auto theta = acos(-p.y);
        auto phi = atan2(-p.z, p.x) + Pi;

        u = phi / (2 * Pi);
        v = theta / Pi;
    }
};

inline bool sphere::hit(const ray& r, Float t_min, Float t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center;
    auto a = r.direction().LengthSquared();
    auto half_b = Dot(oc, r.direction());
    auto c = oc.LengthSquared() - radius * radius;

    auto discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return false;
    auto sqrtd = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    auto root = (-half_b - sqrtd) / a;
    if (root < t_min || t_max < root) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || t_max < root)
            return false;
    }

    rec.time = root;
    rec.p = r.at(rec.time);
    vec3 outward_normal = (rec.p - center) / radius;
    rec.set_face_normal(r, outward_normal);
    point3 p = point3(outward_normal.x, outward_normal.y, outward_normal.z);
    get_sphere_uv(p, rec.u, rec.v);
    rec.mat_ptr = mat_ptr;

    return true;
}

inline bool sphere::bounding_box(Float time0, Float time1, aabb& output_box) const {
    output_box = aabb(
        center - vec3(radius, radius, radius),
        center + vec3(radius, radius, radius));
    return true;
}

class Sphere : public Shape {
public:
    // Sphere Public Methods
    /*zmin zmax给出截断范围，均为局部坐标。phi为[0,2*Pi]*/
    Sphere(shared_ptr<Transform> ObjectToWorld, shared_ptr<Transform> WorldToObject,
        bool reverseOrientation, Float radius, Float zMin, Float zMax,
        Float phiMax)
        : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
        radius(radius),
        zMin(Clamp(std::min(zMin, zMax), -radius, radius)),
        zMax(Clamp(std::max(zMin, zMax), -radius, radius)),
        thetaMin(std::acos(Clamp(std::min(zMin, zMax) / radius, -1, 1))),
        thetaMax(std::acos(Clamp(std::max(zMin, zMax) / radius, -1, 1))),
        phiMax(Radians(Clamp(phiMax, 0, 360))) {}

    Sphere(shared_ptr<Transform> ObjectToWorld, shared_ptr<Transform> WorldToObject, Float radius)
        :Shape(ObjectToWorld,WorldToObject,false),
        radius(radius),
        zMin(-radius),zMax(radius),thetaMin(0),thetaMax(Pi),phiMax(2*Pi){}

    Bounds3f ObjectBound() const;
    bool Intersect(const Ray& ray, Float* tHit, SurfaceInteraction* isect,
        bool testAlphaTexture) const;
    bool IntersectP(const Ray& ray, bool testAlphaTexture) const;
    Float Area() const;
    //Interaction Sample(const Point2f& u, Float* pdf) const;
    //Interaction Sample(const Interaction& ref, const Point2f& u,
    //    Float* pdf) const;
    //Float Pdf(const Interaction& ref, const Vector3f& wi) const;
    //Float SolidAngle(const Point3f& p, int nSamples) const;

private:
    // Sphere Private Data
    const Float radius;
    const Float zMin, zMax;
    const Float thetaMin, thetaMax, phiMax;
};

//std::shared_ptr<Shape> CreateSphereShape(const Transform* o2w,
//    const Transform* w2o,
//    bool reverseOrientation,
//    const ParamSet& params);

#endif