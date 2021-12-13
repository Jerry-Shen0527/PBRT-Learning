#ifndef AARECT_H
#define AARECT_H

#include "rtweekend.h"

#include "primitive.h"

class xy_rect : public Primitive {
public:
    xy_rect() {}

    xy_rect(double _x0, double _x1, double _y0, double _y1, double _k,
        shared_ptr<material> mat)
        : x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mp(mat) {};

    virtual bool Intersect(const ray& r, double t_min, double t_max, SurfaceInteraction& inter_record) const override;
    virtual bool IntersectP(const ray& r, double t_min, double t_max) const override;
    virtual aabb bounding_box() const override {
        // The bounding box must have non-zero width in each dimension, so pad the Z
        // dimension a small amount.
        
        return aabb(Point3f(x0, y0, k - 0.0001), Point3f(x1, y1, k + 0.0001));
        //return true;
    }

public:
    shared_ptr<material> mp;
    double x0, x1, y0, y1, k;
};

class xz_rect : public Primitive {
public:
    xz_rect() {}

    xz_rect(double _x0, double _x1, double _z0, double _z1, double _k,
        shared_ptr<material> mat)
        : x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(mat) {};

    virtual bool Intersect(const ray& r, double t_min, double t_max, SurfaceInteraction& inter_record) const override;
    virtual bool IntersectP(const ray& r, double t_min, double t_max) const override;
    virtual aabb bounding_box() const override {
        // The bounding box must have non-zero width in each dimension, so pad the Y
        // dimension a small amount.
        return aabb(Point3f(x0, k - 0.0001, z0), Point3f(x1, k + 0.0001, z1));
        //return aabb(Point3f(x0, k, z0), Point3f(x1, k, z1));
    }

    virtual double pdf_value(const Point3f& origin, const Vector3f& v) const {
        SurfaceInteraction rec;
        if (!this->Intersect(ray(origin, v), 0.001, infinity, rec))
            return 0;

        auto area = (x1 - x0) * (z1 - z0);
        //cout << rec.t << "1" << v.LengthSquared() << endl;
        auto distance_squared = rec.t * rec.t * v.LengthSquared();
        //std::cout << "2 " << rec.t << std::endl;
        auto cosine = fabs(Dot(v, rec.normal) / v.Length());
        //std::cout << distance_squared<<" "<<(cosine)<<" "<<  area << std::endl;
        return distance_squared / (cosine * area);
    }

    virtual Vector3f random(const Point3f& origin) const {
        auto random_point = Point3f(random_double(x0, x1), k, random_double(z0, z1));
        return Vector3f(random_point - origin);
    }

public:
    shared_ptr<material> mp;
    double x0, x1, z0, z1, k;
};

class yz_rect : public Primitive {
public:
    yz_rect() {}

    yz_rect(double _y0, double _y1, double _z0, double _z1, double _k,
        shared_ptr<material> mat)
        : y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat) {};

    virtual bool Intersect(const ray& r, double t_min, double t_max, SurfaceInteraction& inter_record) const override;
    virtual bool IntersectP(const ray& r, double t_min, double t_max) const override;
    virtual aabb bounding_box() const override {
        // The bounding box must have non-zero width in each dimension, so pad the X
        // dimension a small amount.
        return aabb(Point3f(k - 0.0001, y0, z0), Point3f(k + 0.0001, y1, z1));
        //return true;
    }

public:
    shared_ptr<material> mp;
    double y0, y1, z0, z1, k;
};

#endif