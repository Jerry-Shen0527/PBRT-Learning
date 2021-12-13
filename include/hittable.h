#ifndef HITTABLE_H
#define HITTABLE_H

#include "ray.h"

#include "rtweekend.h"
#include "aabb.h"

class material;

struct hit_record {
    Point3f p;
    Normal3f normal;
    double t;
    shared_ptr<material> mat_ptr;
    double u;
    double v;
    bool front_face;

    inline void set_face_normal(const ray& r, const Normal3f& outward_normal) {
        front_face = Dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

class hittable {
public:
    //pure virtual function, programmers who specification inherits this class must implement this function.
    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
    virtual bool bounding_box(double time0, double time1, aabb& output_box) const = 0;

    virtual double pdf_value(const Point3f& o, const Vector3f& v) const {
        return 0.0;
    }

    virtual Vector3f random(const Point3f& o) const {
        return Vector3f(1, 0, 0);
    }
};

class translate : public hittable {
public:
    translate(shared_ptr<hittable> p, const Vector3f& displacement)
        : ptr(p), offset(displacement) {}

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override;

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override;

public:
    shared_ptr<hittable> ptr;
    Vector3f offset;
};

class rotate_y : public hittable {
public:
    rotate_y(shared_ptr<hittable> p, double angle);

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override;

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        output_box = bbox;
        return hasbox;
    }

public:
    shared_ptr<hittable> ptr;
    double sin_theta;
    double cos_theta;
    bool hasbox;
    aabb bbox;
};

class flip_face : public hittable {
public:
    flip_face(shared_ptr<hittable> p) : ptr(p) {}

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override {

        if (!ptr->hit(r, t_min, t_max, rec))
            return false;

        rec.front_face = !rec.front_face;
        return true;
    }

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        return ptr->bounding_box(time0, time1, output_box);
    }

public:
    shared_ptr<hittable> ptr;
};

#endif
