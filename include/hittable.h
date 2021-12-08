#ifndef HITTABLE_H
#define HITTABLE_H

#include "rtweekend.h"
#include "aabb.h"
class material;
class Shape;
class aabb;
class RGBSpectrum;
class Primitive;

struct hit_record {
    point3 p;
    vec3 normal;
    Float time;
    Float u;
    Float v;
    bool front_face;
    vec3 pError;
    Normal n;
    vec3 wo;
    shared_ptr<material> mat_ptr;
    hit_record():time(0){}
    hit_record(const point3& p, const Normal& n, const vec3& pError,
        const vec3& wo, Float time)
        : p(p), time(time), pError(pError), wo(wo), n(n){ }

    inline void set_face_normal(const ray& r, const vec3& outward_normal) {
        front_face = Dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
    bool IsSurfaceInteraction() const {
        return n != Normal();
    }
};

typedef hit_record Interaction;

class SurfaceInteraction :public hit_record {
public:

    struct {
        Normal n;
        vec3 dpdu, dpdv;
        Normal dndu, dndv;
    } shading;

    SurfaceInteraction() {}
    SurfaceInteraction(
        const point3& p, const vec3& pError, Point2f uv,
        const vec3& wo, const vec3& dpdu, const vec3& dpdv,
        const Normal& dndu, const Normal& dndv, Float time, const Shape* shape,
        int faceIndex = 0);

    void SetShadingGeometry(const vec3& dpdus,
        const vec3& dpdvs, const Normal& dndus,
        const Normal& dndvs, bool orientationIsAuthoritative);

    RGBSpectrum Le(const Vector3f& w) const;

    Point2f uv;
    vec3 dpdu, dpdv;
    Normal dndu, dndv;
    const Shape* shape = nullptr;
    const Primitive* primitive = nullptr;
    //BSDF* bsdf = nullptr;
    //BSSRDF* bssrdf = nullptr;

    mutable vec3 dpdx, dpdy;
    mutable Float dudx = 0, dvdx = 0, dudy = 0, dvdy = 0;

};


class hittable {
public:
    virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec) const = 0;

    virtual bool bounding_box(Float time0, Float time1, aabb& output_box) const = 0;
    virtual Float pdf_value(const point3& o, const vec3& v) const {
        return 0.0;
    }

    virtual vec3 random(const vec3& o) const {
        return vec3(1, 0, 0);
    }
};

class translate : public hittable {
public:
    translate(shared_ptr<hittable> p, const vec3& displacement)
        : ptr(p), offset(displacement) {}

    virtual bool hit(
        const ray& r, Float t_min, Float t_max, hit_record& rec) const override;

    virtual bool bounding_box(Float time0, Float time1, aabb& output_box) const override;

public:
    shared_ptr<hittable> ptr;
    vec3 offset;
};

class rotate_y : public hittable {
public:
    rotate_y(shared_ptr<hittable> p, Float angle);

    virtual bool hit(
        const ray& r, Float t_min, Float t_max, hit_record& rec) const override;

    virtual bool bounding_box(Float time0, Float time1, aabb& output_box) const override;

public:
    shared_ptr<hittable> ptr;
    Float sin_theta;
    Float cos_theta;
    bool hasbox;
    aabb bbox;
};



#endif