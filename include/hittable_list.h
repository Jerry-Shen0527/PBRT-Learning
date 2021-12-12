#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H

#include "hittable.h"
#include "rtweekend.h"
#include "primitive.h"

#include <memory>
#include <vector>
#include <algorithm>

using std::shared_ptr;
using std::make_shared;

class hittable_list : public hittable {
public:
    hittable_list() {}
    hittable_list(shared_ptr<hittable> object) { add(object); }

    void clear() { objects.clear(); }
    void add(shared_ptr<hittable> object) { objects.push_back(object); }
    void setTime(Float t0, Float t1) { _time0 = t0; _time1 = t1; };

    virtual bool hit(
        const ray& r, Float t_min, Float t_max ,hit_record& rec) const override;

    virtual bool bounding_box(
        Float time0, Float time1, aabb& output_box) const override;
public:
    std::vector<shared_ptr<hittable>> objects;
    Float _time0, _time1;
};

bool hittable_list::hit(const ray& r, Float t_min, Float t_max, hit_record& rec) const {
    hit_record temp_rec;
    bool hit_anything = false;
    auto closest_so_far = t_max;

    for (const auto& object : objects) {
        if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.time;
            rec = temp_rec;
        }
    }
    return hit_anything;
}

//list中所有元素的bounding box，若有元素无bounding box 返回false
bool hittable_list::bounding_box(Float time0, Float time1, aabb& output_box) const {
    if (objects.empty()) return false;

    aabb temp_box;
    bool first_box = true;

    for (const auto& object : objects) {
        if (!object->bounding_box(time0, time1, temp_box)) return false;
        output_box = first_box ? temp_box : Union(output_box, temp_box);
        first_box = false;
    }

    return true;
}
#endif

#ifndef BVH_H
#define BVH_H

class bvh_node : public hittable {
public:
    bvh_node();

    bvh_node(const hittable_list& list, Float time0, Float time1)
        : bvh_node(list.objects, 0, list.objects.size(), time0, time1)
    {}

    bvh_node(
        const std::vector<shared_ptr<hittable>>& src_objects,
        size_t start, size_t end, Float time0, Float time1);

    virtual bool hit(
        const ray& r, Float t_min, Float t_max, hit_record& rec) const override;

    virtual bool bounding_box(Float time0, Float time1, aabb& output_box) const override;

public:
    shared_ptr<hittable> left;
    shared_ptr<hittable> right;
    aabb box;
};

inline bool box_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b, int axis) {
    aabb box_a;
    aabb box_b;

    if (!a->bounding_box(0, 0, box_a) || !b->bounding_box(0, 0, box_b))
        std::cerr << "No bounding box in bvh_node constructor.\n";

    return box_a.min()[axis] < box_b.min()[axis];
}


bool box_x_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return box_compare(a, b, 0);
}

bool box_y_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return box_compare(a, b, 1);
}

bool box_z_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return box_compare(a, b, 2);
}



bvh_node::bvh_node(
    const std::vector<shared_ptr<hittable>>& src_objects,
    size_t start, size_t end, Float time0, Float time1
) {
    auto objects = src_objects; // Create a modifiable array of the source scene objects

    int axis = RandomInt(0, 2);
    auto comparator = (axis == 0) ? box_x_compare
        : (axis == 1) ? box_y_compare
        : box_z_compare;

    size_t object_span = end - start;

    if (object_span == 1) {
        left = right = objects[start];
    }
    else if (object_span == 2) {
        if (comparator(objects[start], objects[start + 1])) {
            left = objects[start];
            right = objects[start + 1];
        }
        else {
            left = objects[start + 1];
            right = objects[start];
        }
    }
    else {
        std::sort(objects.begin() + start, objects.begin() + end, comparator);

        auto mid = start + object_span / 2;
        left = make_shared<bvh_node>(objects, start, mid, time0, time1);
        right = make_shared<bvh_node>(objects, mid, end, time0, time1);
    }

    aabb box_left, box_right;

    if (!left->bounding_box(time0, time1, box_left)
        || !right->bounding_box(time0, time1, box_right)
        )
        std::cerr << "No bounding box in bvh_node constructor.\n";

    box = Union(box_left, box_right);
}

bool bvh_node::bounding_box(Float time0, Float time1, aabb& output_box) const {
    output_box = box;
    return true;
}

bool bvh_node::hit(const ray& r, Float t_min, Float t_max, hit_record& rec) const {
    if (!box.hit(r, t_min, t_max))
        return false;

    bool hit_left = left->hit(r, t_min, t_max, rec);
    bool hit_right = right->hit(r, t_min, hit_left ? rec.time : t_max, rec);

    return hit_left || hit_right;
}


inline bool box_compare(const shared_ptr<Primitive> a, const shared_ptr<Primitive> b, int axis) {
    aabb box_a;
    aabb box_b;
    return a->WorldBound().min()[axis] < b->WorldBound().max()[axis];
}

class pbrt_bvh_node :public Primitive
{
public:
    pbrt_bvh_node(const std::vector<shared_ptr<Primitive>> objects,int start,int end);
    ~pbrt_bvh_node();

private:

    // 通过 Primitive 继承
    virtual aabb WorldBound() const override;

    virtual bool Intersect(const Ray& r, SurfaceInteraction* isect) const override;

    virtual bool IntersectP(const Ray& r) const override;

    virtual const Material* GetMaterial() const override;

    virtual void ComputeScatteringFunctions(SurfaceInteraction* isect, TransportMode mode, bool allowMultipleLobes) const override;

    shared_ptr<Primitive> left;
    shared_ptr<Primitive> right;
    aabb box;
};

pbrt_bvh_node::pbrt_bvh_node(const std::vector<shared_ptr<Primitive>> objects, int start, int end)
{
    //auto objects = src_objects; // Create a modifiable array of the source scene objects

    //int axis = RandomInt(0, 2);
    //auto comparator = (axis == 0) ? box_x_compare
    //    : (axis == 1) ? box_y_compare
    //    : box_z_compare;

    //size_t object_span = end - start;
    int object_span = end - start;
    if (object_span == 1) {
        left = right = objects[start];
    }
    //else if (object_span == 2) {
    //    if (comparator(objects[start], objects[start + 1])) {
    //        left = objects[start];
    //        right = objects[start + 1];
    //    }
    //    else {
    //        left = objects[start + 1];
    //        right = objects[start];
    //    }
    //}
    //else {
    //    std::sort(objects.begin() + start, objects.begin() + end, comparator);

    //    auto mid = start + object_span / 2;
    //    left = make_shared<bvh_node>(objects, start, mid, time0, time1);
    //    right = make_shared<bvh_node>(objects, mid, end, time0, time1);
    //}

    //aabb box_left, box_right;

    //if (!left->bounding_box(time0, time1, box_left)
    //    || !right->bounding_box(time0, time1, box_right)
    //    )
    //    std::cerr << "No bounding box in bvh_node constructor.\n";

    //box = Union(box_left, box_right);


}

pbrt_bvh_node::~pbrt_bvh_node()
{
}

aabb pbrt_bvh_node::WorldBound() const
{
    return box;
}

bool pbrt_bvh_node::Intersect(const Ray& r, SurfaceInteraction* isect) const
{

    if (!box.hit(r, r.tMin, r.tMax))
           return false;

    bool hit_left = left->Intersect(r, isect);
    bool hit_right = right->Intersect(r, isect);

    return hit_left || hit_right;
}

bool pbrt_bvh_node::IntersectP(const Ray& r) const
{
    if (!box.hit(r, r.tMin, r.tMax))
        return false;

    bool hit_left = left->IntersectP(r);
    bool hit_right = right->IntersectP(r);

    return hit_left || hit_right;
}

const Material* pbrt_bvh_node::GetMaterial() const
{
    return nullptr;
}

void pbrt_bvh_node::ComputeScatteringFunctions(SurfaceInteraction* isect, TransportMode mode, bool allowMultipleLobes) const
{
}

#endif