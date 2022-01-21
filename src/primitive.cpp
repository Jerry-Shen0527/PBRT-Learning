


// core/primitive.cpp*
#include "primitive.h"
//#include "light.h"
#include "interaction.h"
//#include "stats.h"

namespace pbrt {

//STAT_MEMORY_COUNTER("Memory/Primitives", primitiveMemory);

// Primitive Method Definitions
Primitive::~Primitive() {}

//Aggregate
/*
const AreaLight *Aggregate::GetAreaLight() const {
    //LOG(FATAL) <<
      //  "Aggregate::GetAreaLight() method"
        //"called; should have gone to GeometricPrimitive";
    return nullptr;
}

const material *Aggregate::GetMaterial() const {
    LOG(FATAL) <<
        "Aggregate::GetMaterial() method"
        "called; should have gone to GeometricPrimitive";
    return nullptr;
}

void Aggregate::ComputeScatteringFunctions(SurfaceInteraction *isect,
                                           MemoryArena &arena,
                                           TransportMode mode,
                                           bool allowMultipleLobes) const {
    LOG(FATAL) <<
        "Aggregate::ComputeScatteringFunctions() method"
        "called; should have gone to GeometricPrimitive";
}*/

// TransformedPrimitive Method Definitions
TransformedPrimitive::TransformedPrimitive(std::shared_ptr<Primitive> &primitive,
                                           const AnimatedTransform &PrimitiveToWorld)
    : primitive(primitive), PrimitiveToWorld(PrimitiveToWorld) {
    //primitiveMemory += sizeof(*this);
}

bool TransformedPrimitive::Intersect(const Ray &r,
                                     SurfaceInteraction *isect) const {
    // Compute _ray_ after transformation by _PrimitiveToWorld_
    Transform InterpolatedPrimToWorld;
    PrimitiveToWorld.Interpolate(r.time, &InterpolatedPrimToWorld);
    Ray ray = Inverse(InterpolatedPrimToWorld)(r);
    if (!primitive->Intersect(ray, isect)) return false;
    r.tMax = ray.tMax;
    // Transform instance's intersection data to world space
    if (!InterpolatedPrimToWorld.IsIdentity())
        *isect = InterpolatedPrimToWorld(*isect);
    CHECK_GE(Dot(isect->n, isect->shading.n), 0);
    return true;
}

bool TransformedPrimitive::IntersectP(const Ray &r) const {
    Transform InterpolatedPrimToWorld;
    PrimitiveToWorld.Interpolate(r.time, &InterpolatedPrimToWorld);
    Transform InterpolatedWorldToPrim = Inverse(InterpolatedPrimToWorld);
    return primitive->IntersectP(InterpolatedWorldToPrim(r));
}

// GeometricPrimitive Method Definitions
GeometricPrimitive::GeometricPrimitive(const std::shared_ptr<Shape> &shape,
                                       const std::shared_ptr<Material> &material,
                                       const std::shared_ptr<AreaLight> &areaLight)
    : shape(shape),
    material_(material),
    areaLight(areaLight){
    //std::cout << "GeometricPrimitive::GeometricPrimitive material_: " << material_ << std::endl;
    //primitiveMemory += sizeof(*this);
}

Bounds3f GeometricPrimitive::WorldBound() const { return shape->WorldBound(); }

bool GeometricPrimitive::IntersectP(const Ray &r) const {
    return shape->IntersectP(r);
}

bool GeometricPrimitive::Intersect(const Ray &r,
                                   SurfaceInteraction *isect) const {
    Float tHit;
    if (!shape->Intersect(r, &tHit, isect)) return false;
    r.tMax = tHit;
    isect->primitive = this;
    //add material
    //std::cout << "GeometricPrimitive::Intersect material_: " << material_ << std::endl;
    isect->mat_ptr = material_;
    CHECK_GE(Dot(isect->n, isect->shading.n), 0.);
    // Initialize _SurfaceInteraction::mediumInterface_ after _Shape_
    // intersection    
    //if (mediumInterface.IsMediumTransition())
      //  isect->mediumInterface = mediumInterface;
    //else
      //  isect->mediumInterface = MediumInterface(r.medium);
    return true;
}

const AreaLight *GeometricPrimitive::GetAreaLight() const {
    return areaLight.get();
}

const Material *GeometricPrimitive::GetMaterial() const {
    return material_.get();
}

void GeometricPrimitive::ComputeScatteringFunctions(
    SurfaceInteraction *isect, MemoryArena &arena, TransportMode mode,
    bool allowMultipleLobes) const {
    //ProfilePhase p(Prof::ComputeScatteringFuncs);
    if (material_)
        material_->ComputeScatteringFunctions(isect, arena, mode,
                                             allowMultipleLobes);
    CHECK_GE(Dot(isect->n, isect->shading.n), 0.);
}

Bounds3f PrimitiveLists::WorldBound() const {
    Bounds3f worBound;
    for (auto& primitive : primitives) {
        worBound = Union(worBound, primitive->WorldBound());
    }
    return worBound;
}
bool PrimitiveLists::Intersect(const Ray& r, SurfaceInteraction* isect) const {
    bool hit_prims = false;
    for (auto& primitive : primitives) {
        if (primitive->IntersectP(r))
        {
            hit_prims = true;
            primitive->Intersect(r, isect);
        }
    }
    return hit_prims;
}
bool PrimitiveLists::IntersectP(const Ray& r) const {
    for (auto& primitive : primitives) {
        if (primitive->IntersectP(r))
            return true;
    }
    return false;
}


}  // namespace pbrt
