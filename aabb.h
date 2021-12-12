#ifndef AABB_H
#define AABB_H

#include "rtweekend.h"
#include "mtransform.h"
#include "assert.h"
//我修改过的
class aabb {
public:
    aabb() {}
    aabb(Transform otw, Transform wto, const point3& min_, const point3& max_) {
        minimum = min_; maximum = max_;
        ObjectToWorld = otw;
        WorldToObject = wto;
    }

    aabb(const point3& min_, const point3& max_) {
        minimum = min_; maximum = max_;
    }


    point3 min() const { return minimum; }
    point3 max() const { return maximum; }

    bool hit(const ray& r, double t_min, double t_max) const;

    point3 minimum;
    point3 maximum;
    Transform ObjectToWorld, WorldToObject;
};

inline bool aabb::hit(const ray& wr, double t_min, double t_max) const {
    ray r = WorldToObject.ray_transform(wr);
    for (int a = 0; a < 3; a++) {

        auto invD = 1.0f / r.direction()[a];
        auto t0 = (min()[a] - r.origin()[a]) * invD;
        auto t1 = (max()[a] - r.origin()[a]) * invD;
        if (invD < 0.0f)
            std::swap(t0, t1);
        t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;
        if (t_max <= t_min)
            return false;
    }
    return true;
}

aabb surrounding_box(aabb box0, aabb box1) {

    assert(box0.ObjectToWorld == box1.ObjectToWorld);

    point3 small(fmin(box0.min().x(), box1.min().x()),
        fmin(box0.min().y(), box1.min().y()),
        fmin(box0.min().z(), box1.min().z()));

    point3 big(fmax(box0.max().x(), box1.max().x()),
        fmax(box0.max().y(), box1.max().y()),
        fmax(box0.max().z(), box1.max().z()));

    return aabb(box0.ObjectToWorld,box0.WorldToObject,small, big);
}

#endif
