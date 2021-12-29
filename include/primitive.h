
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_PRIMITIVE_H
#define PBRT_CORE_PRIMITIVE_H

// core/primitive.h*
//#include "pbrt.h"
#include "shape.h"
#include "material.h"
#include "medium.h"
#include "transform.h"


//delete MediumInterface and ComputeScatteringFunctions, modify Material-->material

namespace pbrt {
    class AreaLight;
// Primitive Declarations
class Primitive {
  public:
    // Primitive Interface
    virtual ~Primitive();
    virtual Bounds3f WorldBound() const = 0;
    virtual bool Intersect(const Ray &r, SurfaceInteraction *) const = 0;
    virtual bool IntersectP(const Ray &r) const = 0;
    virtual const AreaLight *GetAreaLight() const = 0;
    virtual const Material *GetMaterial() const = 0;
    //virtual const material* GetMaterial() const = 0;
    
    virtual void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                            MemoryArena &arena,
                                            TransportMode mode,
                                            bool allowMultipleLobes) const = 0;
};

// GeometricPrimitive Declarations
class GeometricPrimitive : public Primitive {
  public:
    // GeometricPrimitive Public Methods
    virtual Bounds3f WorldBound() const;
    virtual bool Intersect(const Ray &r, SurfaceInteraction *isect) const;
    virtual bool IntersectP(const Ray &r) const;
    /*
    GeometricPrimitive(const std::shared_ptr<Shape> &shape,
                       const std::shared_ptr<Material> &material,
                       const std::shared_ptr<AreaLight> &areaLight,
                       const MediumInterface &mediumInterface);*/
    GeometricPrimitive(const std::shared_ptr<Shape>& shape,
        const std::shared_ptr<Material>& material,
        const std::shared_ptr<AreaLight>& areaLight);
    const AreaLight *GetAreaLight() const;
    const Material *GetMaterial() const;
    
    void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                    MemoryArena &arena, TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // GeometricPrimitive Private Data
    std::shared_ptr<Shape> shape;
    //std::shared_ptr<material> material_;
    std::shared_ptr<Material> material_;
    std::shared_ptr<AreaLight> areaLight;
    //MediumInterface mediumInterface;
};

//geo default (enum)
class PrimitiveLists :public Primitive {
public:
    PrimitiveLists(){}
    void add(const std::shared_ptr<Shape>& shape,
        const std::shared_ptr<Material>& material,
        const std::shared_ptr<AreaLight>& areaLight) {
        primitives.push_back(make_shared<GeometricPrimitive>(shape, material, areaLight));
    }
    void add(const std::vector<std::shared_ptr<Shape>>& shapes,
        const std::shared_ptr<Material>& material,
        const std::shared_ptr<AreaLight>& areaLight) {
        for (auto& shape : shapes)
            primitives.push_back(make_shared<GeometricPrimitive>(shape, material, areaLight));
    }
    void add(const std::shared_ptr<Primitive>& prim) {
        primitives.push_back(prim);
    }
    void add(const std::vector<std::shared_ptr<Primitive>>& prims) {
        for (auto& prim : prims)
            primitives.push_back(prim);
    }

    Bounds3f WorldBound() const;
    bool Intersect(const Ray& r, SurfaceInteraction* isect) const;
    bool IntersectP(const Ray& r) const;
    const AreaLight* GetAreaLight() const { return nullptr; }
    const Material* GetMaterial() const { return nullptr; }

    void ComputeScatteringFunctions(SurfaceInteraction* isect,
        MemoryArena& arena, TransportMode mode,
        bool allowMultipleLobes) const {
        std::cerr <<
            "PrimitiveLists::ComputeScatteringFunctions() shouldn't be "
            "called";
    }

public:
    std::vector<shared_ptr<Primitive>> primitives;
};


// TransformedPrimitive Declarations
class TransformedPrimitive : public Primitive {
  public:
    // TransformedPrimitive Public Methods
    TransformedPrimitive(std::shared_ptr<Primitive> &primitive,
                         const AnimatedTransform &PrimitiveToWorld);
    bool Intersect(const Ray &r, SurfaceInteraction *in) const;
    bool IntersectP(const Ray &r) const;
    const AreaLight *GetAreaLight() const { return nullptr; }
    const Material *GetMaterial() const { return nullptr; }
    
    void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                    MemoryArena &arena, TransportMode mode,
                                    bool allowMultipleLobes) const {
        std::cerr <<
            "TransformedPrimitive::ComputeScatteringFunctions() shouldn't be "
            "called";
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
    const AreaLight *GetAreaLight() const;
    const Material *GetMaterial() const;
    /*
    void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                    MemoryArena &arena, TransportMode mode,
                                    bool allowMultipleLobes) const;*/
};

}  // namespace pbrt

#endif  // PBRT_CORE_PRIMITIVE_H
