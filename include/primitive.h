#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include "ray.h"
#include "material.h"
#include "transform.h"

class aabb;
class SurfaceInteraction;
class AreaLight;


class Primitive {
public:
    // Primitive Interface
    virtual ~Primitive(){}
    virtual aabb WorldBound() const = 0;
    virtual bool Intersect(const Ray& r, SurfaceInteraction* isect) const = 0;
    virtual bool IntersectP(const Ray& r) const = 0;
    //virtual const AreaLight* GetAreaLight() const = 0;
    virtual const Material* GetMaterial() const = 0;
    virtual void ComputeScatteringFunctions(SurfaceInteraction* isect,
        //MemoryArena& arena,
        TransportMode mode,
        bool allowMultipleLobes) const = 0;
};

// GeometricPrimitive Declarations
class GeometricPrimitive : public Primitive {
public:
    // GeometricPrimitive Public Methods
    virtual aabb WorldBound() const;
    virtual bool Intersect(const Ray& r, SurfaceInteraction* isect) const;
    virtual bool IntersectP(const Ray& r) const;
    GeometricPrimitive(const std::shared_ptr<Shape>& shape,
        const std::shared_ptr<Material>& material,
        const std::shared_ptr<AreaLight>& areaLight//,
        //const MediumInterface& mediumInterface
    );
    GeometricPrimitive(const std::shared_ptr<Shape>& shape);
    const AreaLight* GetAreaLight() const;
    const Material* GetMaterial() const;
    void ComputeScatteringFunctions(SurfaceInteraction* isect,
        //MemoryArena& arena, 
        TransportMode mode,
        bool allowMultipleLobes) const;

private:
    // GeometricPrimitive Private Data
    std::shared_ptr<Shape> shape;
    std::shared_ptr<Material> material;
    std::shared_ptr<AreaLight> areaLight;
    //MediumInterface mediumInterface;
};

// TransformedPrimitive Declarations
class TransformedPrimitive : public Primitive {
public:
    // TransformedPrimitive Public Methods
    TransformedPrimitive(std::shared_ptr<Primitive>& primitive,
        const AnimatedTransform& PrimitiveToWorld);
    bool Intersect(const Ray& r, SurfaceInteraction* in) const;
    bool IntersectP(const Ray& r) const;
    const AreaLight* GetAreaLight() const { return nullptr; }
    const Material* GetMaterial() const { return nullptr; }
    void ComputeScatteringFunctions(SurfaceInteraction* isect,
        /*MemoryArena& arena,*/ TransportMode mode,
        bool allowMultipleLobes) const {
        //LOG(FATAL) <<
        //    "TransformedPrimitive::ComputeScatteringFunctions() shouldn't be "
        //    "called";
    }
    Bounds3f WorldBound() const {
        return PrimitiveToWorld.MotionBounds(primitive->WorldBound());
    }

private:
    // TransformedPrimitive Private Data
    std::shared_ptr<Primitive> primitive;
    const AnimatedTransform PrimitiveToWorld;
};


// Aggregate Declarations
class Aggregate : public Primitive {
public:
    // Aggregate Public Methods
    const AreaLight* GetAreaLight() const;
    const Material* GetMaterial() const;
    void ComputeScatteringFunctions(SurfaceInteraction* isect,
        /*MemoryArena& arena,*/ TransportMode mode,
        bool allowMultipleLobes) const;
};


#endif // !PRIMITIVE_H
