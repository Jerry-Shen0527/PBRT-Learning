#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif


#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
#include "vec3.h"
#include "pdf.h"

class sphere : public hittable {
public:
    sphere() {}
    sphere(Point3f cen, Float r, shared_ptr<material> m)
        : center(cen), radius(r), mat_ptr(m) {};

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override;
    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override;
    virtual double pdf_value(const Point3f& o, const Vector3f& v) const override;
    virtual Vector3f random(const Point3f& o) const override;

private:
    static void get_sphere_uv(const Vector3f& p, double& u, double& v);

public:
    Point3f center;
    double radius;
    shared_ptr<material> mat_ptr;
};

#include "shape.h"

namespace pbrt {

    // Sphere Declarations
    class Sphere : public Shape {
    public:
        // Sphere Public Methods
        Sphere(const Transform* ObjectToWorld, const Transform* WorldToObject,
            bool reverseOrientation, Float radius, Float zMin, Float zMax,
            Float phiMax)
            : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
            radius(radius),
            zMin(Clamp(std::min(zMin, zMax), -radius, radius)),
            zMax(Clamp(std::max(zMin, zMax), -radius, radius)),
            thetaMin(std::acos(Clamp(std::min(zMin, zMax) / radius, -1, 1))),
            thetaMax(std::acos(Clamp(std::max(zMin, zMax) / radius, -1, 1))),
            phiMax(Radians(Clamp(phiMax, 0, 360))) {}
        Bounds3f ObjectBound() const;
        bool Intersect(const Ray& ray, Float* tHit, SurfaceInteraction* isect,
            bool testAlphaTexture) const;
        bool IntersectP(const Ray& ray, bool testAlphaTexture) const;
        Float Area() const;
        Interaction Sample(const Point2f& u, Float* pdf) const;
        Interaction Sample(const Interaction& ref, const Point2f& u,
            Float* pdf) const;
        Float Pdf(const Interaction& ref, const Vector3f& wi) const;
        Float SolidAngle(const Point3f& p, int nSamples) const;

    private:
        // Sphere Private Data
        const Float radius;
        const Float zMin, zMax;
        const Float thetaMin, thetaMax, phiMax;
    };
    /*
    std::shared_ptr<Shape> CreateSphereShape(const Transform* o2w,
        const Transform* w2o,
        bool reverseOrientation,
        const ParamSet& params);
    */
}  // namespace pbrt



#endif
