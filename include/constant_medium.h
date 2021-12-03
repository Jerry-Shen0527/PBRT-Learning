#ifndef CONSTANT_MEDIUM_H
#define CONSTANT_MEDIUM_H

#include "rtweekend.h"

#include "hittable.h"
#include "material.h"
#include "texture.h"

class constant_medium : public hittable {
public:
    constant_medium(shared_ptr<hittable> b, Float d, shared_ptr<texture> a)
        : boundary(b),
        neg_inv_density(-1 / d),
        phase_function(make_shared<isotropic>(a))
    {}

    constant_medium(shared_ptr<hittable> b, Float d, color c)
        : boundary(b),
        neg_inv_density(-1 / d),
        phase_function(make_shared<isotropic>(c))
    {}

    virtual bool hit(
        const ray& r, Float t_min, Float t_max, hit_record& rec) const override;

    virtual bool bounding_box(Float time0, Float time1, aabb& output_box) const override {
        return boundary->bounding_box(time0, time1, output_box);
    }

public:
    shared_ptr<hittable> boundary;
    shared_ptr<material> phase_function;
    Float neg_inv_density;
};

bool constant_medium::hit(const ray& r, Float t_min, Float t_max, hit_record& rec) const {
    // Print occasional samples when debugging. To enable, set enableDebug true.
    const bool enableDebug = false;
    const bool debugging = enableDebug && RandomFloat() < 0.00001;

    hit_record rec1, rec2;

    if (!boundary->hit(r, -Infinity, Infinity, rec1))
        return false;

    if (!boundary->hit(r, rec1.time + 0.0001, Infinity, rec2))
        return false;

    if (debugging) std::cerr << "\nt_min=" << rec1.time << ", t_max=" << rec2.time << '\n';

    if (rec1.time < t_min) rec1.time = t_min;
    if (rec2.time > t_max) rec2.time = t_max;

    if (rec1.time >= rec2.time)
        return false;

    if (rec1.time < 0)
        rec1.time = 0;

    const auto ray_length = r.direction().Length();
    const auto distance_inside_boundary = (rec2.time - rec1.time) * ray_length;
    const auto hit_distance = neg_inv_density * log(RandomFloat());

    if (hit_distance > distance_inside_boundary)
        return false;

    rec.time = rec1.time + hit_distance / ray_length;
    rec.p = r.at(rec.time);

   /*if (debugging) {
        std::cerr << "hit_distance = " << hit_distance << '\n'
            << "rec.t = " << rec.t << '\n'
            << "rec.p = " << rec.p << '\n';
    }*/ 

    rec.normal = vec3(1, 0, 0);  // arbitrary
    rec.front_face = true;     // also arbitrary
    rec.mat_ptr = phase_function;

    return true;
}
#endif