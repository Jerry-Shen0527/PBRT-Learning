#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H

#include "hittable.h"

#include <memory>
#include <vector>

using std::shared_ptr;
using std::make_shared;

class hittable_list : public hittable {
public:
    hittable_list() {}
    hittable_list(shared_ptr<hittable> object) { add(object); }

    void clear() { objects.clear(); }
    void add(shared_ptr<hittable> object) { objects.push_back(object); }

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override;
    virtual bool bounding_box(
        double time0, double time1, aabb& output_box) const override;

    virtual double pdf_value(const Point3f& o, const Vector3f& v) const override;

    virtual Vector3f random(const Point3f& o) const override;

public:
    std::vector<shared_ptr<hittable>> objects;
};


#endif
