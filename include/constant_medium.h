#ifndef CONSTANT_MEDIUM_H
#define CONSTANT_MEDIUM_H

#include "rtweekend.h"

#include "hittable.h"
#include "material.h"
#include "texture.h"

class constant_medium : public hittable {
public:
    constant_medium(shared_ptr<hittable> b, double d, shared_ptr<texture> a)
        : boundary(b),
        neg_inv_density(-1 / d),
        phase_function(make_shared<isotropic>(a))
    {}

    constant_medium(shared_ptr<hittable> b, double d, color c)
        : boundary(b),
        neg_inv_density(-1 / d),
        phase_function(make_shared<isotropic>(c))
    {}

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override;

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        return boundary->bounding_box(time0, time1, output_box);
    }

public:
    shared_ptr<hittable> boundary;
    shared_ptr<material> phase_function;
    double neg_inv_density;
};

#endif