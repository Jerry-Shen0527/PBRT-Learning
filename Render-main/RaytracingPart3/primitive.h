#pragma once
#include "geo.h"
//#include "material.h"
#include "aabb.h"
#include "ray.h"
using namespace std;
class material;

struct SurfaceInteraction {
    Point3f p;
    Vector3f normal;
    shared_ptr<material> mat_ptr;
    double t;
    bool front_face;
    double u;
    double v;
    Vector3f dpdu, dpdv;
    //Normal3f dndu, dndv;


    inline void set_face_normal(const ray& r, const Vector3f& outward_normal) {
        front_face = Dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

class Primitive
{
public:
    // Primitive Interface
    virtual ~Primitive();
    virtual aabb bounding_box() const = 0;
    virtual bool Intersect(const ray& r, double t_min, double t_max, SurfaceInteraction& inter_record) const = 0;
    virtual bool IntersectP(const ray& r, double t_min, double t_max) const = 0;
    virtual double pdf_value(const Point3f& o, const Vector3f& v) const {
        return 0.0;
    }

    virtual Vector3f random(const Vector3f& o) const {
        return Vector3f(1, 0, 0);
    }
    
    //virtual const material* GetMaterial() const = 0;
};

class translate : public Primitive {
public:
    translate(shared_ptr<Primitive> p, const Vector3f& displacement)
        : ptr(p), offset(displacement) {}

    virtual bool Intersect(const ray& r, double t_min, double t_max, SurfaceInteraction& inter_record) const override;
    virtual bool IntersectP(const ray& r, double t_min, double t_max) const = 0;
    virtual aabb bounding_box() const override;

public:
    shared_ptr<Primitive> ptr;
    Vector3f offset;
};

class rotate_y : public Primitive {
public:
    rotate_y(shared_ptr<Primitive> p, double angle);

    virtual bool Intersect(const ray& r, double t_min, double t_max, SurfaceInteraction& inter_record) const override;
    virtual bool IntersectP(const ray& r, double t_min, double t_max) const = 0;
    virtual aabb bounding_box() const override;

public:
    shared_ptr<Primitive> ptr;
    double sin_theta;
    double cos_theta;
    bool hasbox;
    aabb bbox;
};

class flip_face : public Primitive {
public:
    flip_face(shared_ptr<Primitive> p) : ptr(p) {}

    virtual bool Intersect(const ray& r, double t_min, double t_max, SurfaceInteraction& inter_record) const override {

        if (!ptr->Intersect(r, t_min, t_max, inter_record))
            return false;
        
        inter_record.front_face = !inter_record.front_face;
        return true;
    }
    virtual bool IntersectP(const ray& r, double t_min, double t_max) const override {

        if (!ptr->IntersectP(r, t_min, t_max))
            return false;

        return true;
    }

    virtual aabb bounding_box() const override {
        return ptr->bounding_box();
    }

public:
    shared_ptr<Primitive> ptr;
};

/*class GeometricPrimitive : public Primitive {
public:
    // GeometricPrimitive Public Methods
    virtual bool bounding_box(double time0, double time1, aabb& output_box) const = 0;
    virtual bool Intersect(const ray& r, SurfaceInteraction* isect) const;
    virtual bool IntersectP(const ray& r) const;
    GeometricPrimitive(const std::shared_ptr<Primitive>& primitive,
        const std::shared_ptr<material>& mat);
    
    const material* GetMaterial() const;
    

private:
    // GeometricPrimitive Private Data
    std::shared_ptr<Primitive> primitive;
    std::shared_ptr<material> mat;
    //std::shared_ptr<AreaLight> areaLight;
    //MediumInterface mediumInterface;
};*/

// TransformedPrimitive Declarations
/*class TransformedPrimitive : public Primitive {
public:
    // TransformedPrimitive Public Methods
    TransformedPrimitive(std::shared_ptr<Primitive>& primitive,
        const AnimatedTransform& PrimitiveToWorld);
    bool Intersect(const Ray& r, SurfaceInteraction* in) const;
    bool IntersectP(const Ray& r) const;
    const AreaLight* GetAreaLight() const { return nullptr; }
    const Material* GetMaterial() const { return nullptr; }
    void ComputeScatteringFunctions(SurfaceInteraction* isect,
        MemoryArena& arena, TransportMode mode,
        bool allowMultipleLobes) const {
        LOG(FATAL) <<
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
};*/

// Aggregate Declarations
//class Aggregate : public Primitive {
//public:
//    // Aggregate Public Methods
//    //const AreaLight* GetAreaLight() const;
//    const material* GetMaterial() const;
//    
//};