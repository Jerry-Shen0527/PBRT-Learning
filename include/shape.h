
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_SHAPE_H
#define PBRT_CORE_SHAPE_H

// core/shape.h*
//#include "pbrt.h"
//#include "geometry.h"
#include "interaction.h"
#include "memory.h"
#include "transform.h"

#include "vec3.h"
#include "aabb.h"

namespace pbrt {

// Shape Declarations
class Shape {
  public:
    // Shape Interface
      //default
    Shape(shared_ptr<Transform> ObjectToWorld, shared_ptr<Transform> WorldToObject,
          bool reverseOrientation);
    virtual ~Shape();
    virtual Bounds3f ObjectBound() const = 0;
    virtual Bounds3f WorldBound() const;
    virtual bool Intersect(const Ray &ray, Float *tHit,
                           SurfaceInteraction *isect,
                           bool testAlphaTexture = true) const = 0;
    virtual bool IntersectP(const Ray &ray,
                            bool testAlphaTexture = true) const {
        return Intersect(ray, nullptr, nullptr, testAlphaTexture);
    }
    virtual Float Area() const = 0;
    //Sample
    /*
    // Sample a point on the surface of the shape and return the PDF with
    // respect to area on the surface.
    virtual Interaction Sample(const Point2f &u, Float *pdf) const = 0;
    virtual Float Pdf(const Interaction &) const { return 1 / Area(); }

    // Sample a point on the shape given a reference point |ref| and
    // return the PDF with respect to solid angle from |ref|.
    virtual Interaction Sample(const Interaction &ref, const Point2f &u,
                               Float *pdf) const;
    virtual Float Pdf(const Interaction &ref, const Vector3f &wi) const;

    // Returns the solid angle subtended by the shape w.r.t. the reference
    // point p, given in world space. Some shapes compute this value in
    // closed-form, while the default implementation uses Monte Carlo
    // integration; the nSamples parameter determines how many samples are
    // used in this case.
    virtual Float SolidAngle(const Point3f &p, int nSamples = 512) const;
    */
    // Shape Public Data
    //const Transform *ObjectToWorld, *WorldToObject;
    shared_ptr<Transform> ObjectToWorld, WorldToObject;
    const bool reverseOrientation;
    const bool transformSwapsHandedness;
};

}  // namespace pbrt

#endif  // PBRT_CORE_SHAPE_H
