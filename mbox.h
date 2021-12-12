#ifndef MBOX_H
#define MBOX_H

#include "rtweekend.h"

#include "aarect.h"
#include "hittable_list.h"
#include "mtransform.h"

//box 是重心在原点
class mbox : public hittable {
public:
    mbox() {}
    mbox(Transform otw, Transform wto, const point3& p0, const point3& p1, shared_ptr<material> ptr);

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        output_box = aabb(ObjectToWorld, WorldToObject,box_min, box_max);
        return true;
    }

public:
    point3 box_min;
    point3 box_max;
    hittable_list sides;
    Transform ObjectToWorld, WorldToObject;
};

mbox::mbox(Transform otw, Transform wto, const point3& p0, const point3& p1, shared_ptr<material> ptr) {
    ObjectToWorld = otw;
    WorldToObject = wto;

    box_min = p0;
    box_max = p1;

    sides.add(make_shared<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p1.z(), ptr));
    sides.add(make_shared<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p0.z(), ptr));

    sides.add(make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p1.y(), ptr));
    sides.add(make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p0.y(), ptr));

    sides.add(make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p1.x(), ptr));
    sides.add(make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p0.x(), ptr));
}

bool mbox::hit(const ray& wr, double t_min, double t_max, hit_record& rec) const {
    ray r = WorldToObject.ray_transform(wr);
    bool a= sides.hit(r, t_min, t_max, rec);
    if (!a) return a;
    //r = ObjectToWorld.ray_transform(r);
    rec.p = wr.at(rec.t);
    rec.normal= ObjectToWorld.vector_transform(rec.normal);
    return a;

}

#endif

