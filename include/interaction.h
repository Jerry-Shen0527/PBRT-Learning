


#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_INTERACTION_H
#define PBRT_CORE_INTERACTION_H

// core/interaction.h*
//#include "pbrt.h"
//#include "geometry.h"
#include "vec3.h"
#include "transform.h"
#include "spectrum.h"
#include "material.h"
//#include "medium.h"
//#include "material.h"

//delete medium data  


namespace pbrt {
    //struct MediumInterface;
    class Shape;
    class Primitive;
// Interaction Declarations
struct Interaction {
    // Interaction Public Methods
    Interaction() : time(0) {}
    Interaction(const Point3f &p, const Normal3f &n, const Vector3f &pError,
                const Vector3f &wo, Float time)
        : p(p),
          time(time),
          pError(pError),
          wo(Normalize(wo)),
          n(n){}
    bool IsSurfaceInteraction() const { return n != Normal3f(); }
    ray SpawnRay(const Vector3f &d) const {
        Point3f o = OffsetRayOrigin(p, pError, n, d);
        //return Ray(o, d, infinity, time, GetMedium(d));
        return ray(o, d, time);
    }
    ray SpawnRayTo(const Point3f &p2) const {
        Point3f origin = OffsetRayOrigin(p, pError, n, p2 - p);
        Vector3f d = p2 - p;
        //return ray(origin, d, 1 - ShadowEpsilon, time, GetMedium(d));
        return ray(origin, d, time);
    }
    ray SpawnRayTo(const Interaction &it) const {
        Point3f origin = OffsetRayOrigin(p, pError, n, it.p - p);
        Point3f target = OffsetRayOrigin(it.p, it.pError, it.n, origin - it.p);
        Vector3f d = target - origin;
        //return Ray(origin, d, 1 - ShadowEpsilon, time, GetMedium(d));
        return ray(origin, d, time);
    }
    Interaction(const Point3f &p, const Vector3f &wo, Float time)
        : p(p), time(time), wo(wo) {}
    Interaction(const Point3f &p, Float time)
        : p(p), time(time){}
    bool IsMediumInteraction() const { return !IsSurfaceInteraction(); }
    /*
    const Medium *GetMedium(const Vector3f &w) const {
        return Dot(w, n) > 0 ? mediumInterface.outside : mediumInterface.inside;
    }
    const Medium *GetMedium() const {
        //CHECK_EQ(mediumInterface.inside, mediumInterface.outside);
        return mediumInterface.inside;
    }*/

    // Interaction Public Data
    Point3f p;
    Float time;
    Normal3f n;
    Vector3f pError;
    Vector3f wo;
    //MediumInterface mediumInterface;
};

/*
class MediumInteraction : public Interaction {
  public:
    // MediumInteraction Public Methods
    MediumInteraction() : phase(nullptr) {}
    MediumInteraction(const Point3f &p, const Vector3f &wo, Float time,
                      const Medium *medium, const PhaseFunction *phase)
        : Interaction(p, wo, time, medium), phase(phase) {}
    bool IsValid() const { return phase != nullptr; }

    // MediumInteraction Public Data
    const PhaseFunction *phase;
};
*/

// SurfaceInteraction Declarations
class SurfaceInteraction : public Interaction {
  public:
    // SurfaceInteraction Public Methods
    SurfaceInteraction() {}
    SurfaceInteraction(const Point3f &p, const Vector3f &pError,
                       const Point2f &uv, const Vector3f &wo,
                       const Vector3f &dpdu, const Vector3f &dpdv,
                       const Normal3f &dndu, const Normal3f &dndv, Float time,
                       const Shape *sh,
                       int faceIndex = 0);
    void SetShadingGeometry(const Vector3f &dpdu, const Vector3f &dpdv,
                            const Normal3f &dndu, const Normal3f &dndv,
                            bool orientationIsAuthoritative);
    /*
    void ComputeScatteringFunctions(
        const RayDifferential &ray, MemoryArena &arena,
        bool allowMultipleLobes = false,
        TransportMode mode = TransportMode::Radiance);*/
    void ComputeDifferentials(const RayDifferential &r) const;
    //Spectrum Le(const Vector3f &w) const;

    // SurfaceInteraction Public Data
    Point2f uv;
    Vector3f dpdu, dpdv;
    Normal3f dndu, dndv;
    const Shape *shape = nullptr;
    struct {
        Normal3f n;
        Vector3f dpdu, dpdv;
        Normal3f dndu, dndv;
    } shading;
    const Primitive *primitive = nullptr;
    //BSDF *bsdf = nullptr;
    //BSSRDF *bssrdf = nullptr;
    shared_ptr<material> mat_ptr;
    mutable Vector3f dpdx, dpdy;
    mutable Float dudx = 0, dvdx = 0, dudy = 0, dvdy = 0;

    // Added after book publication. Shapes can optionally provide a face
    // index with an intersection point for use in Ptex texture lookups.
    // If Ptex isn't being used, then this value is ignored.
    int faceIndex = 0;
};

}  // namespace pbrt

#endif  // PBRT_CORE_INTERACTION_H
