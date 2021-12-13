#ifndef MOVING_SPHERE_H
#define MOVING_SPHERE_H

#include "rtweekend.h"

#include "hittable.h"


class moving_sphere : public hittable {
public:
    moving_sphere() {}
    moving_sphere(
        Point3f cen0, Point3f cen1, double _time0, double _time1, double r, shared_ptr<material> m)
        : center0(cen0), center1(cen1), time0(_time0), time1(_time1), radius(r), mat_ptr(m)
    {};

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override;
    virtual bool bounding_box(
        double _time0, double _time1, aabb& output_box) const override;

    Point3f center(double time) const;

public:
    Point3f center0, center1;
    double time0, time1;
    double radius;
    shared_ptr<material> mat_ptr;
};



#endif