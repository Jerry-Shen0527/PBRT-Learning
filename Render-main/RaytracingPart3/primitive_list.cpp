#include "primitive_list.h"
bool Primitive_list::Intersect(const ray& r, double t_min, double t_max, SurfaceInteraction& inter_record) const {
    SurfaceInteraction temp_rec;
    bool hit_anything = false;
    auto closest_so_far = t_max;

    for (const auto& object : objects) {
        if (object->Intersect(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            inter_record = temp_rec;
        }
    }

    return hit_anything;
}

bool Primitive_list::IntersectP(const ray& r, double t_min, double t_max) const {
    SurfaceInteraction temp_rec;
    bool hit_anything = false;
    auto closest_so_far = t_max;

    for (const auto& object : objects) {
        if (object->IntersectP(r, t_min, closest_so_far)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            //inter_record = temp_rec;
        }
    }

    return hit_anything;
}

aabb Primitive_list::bounding_box() const {
    //if (objects.empty()) return nullptr;

    aabb temp_box;
    bool first_box = true;

    for (const auto& object : objects) {
        //if (!object->bounding_box(time0, time1, temp_box)) return false;
        temp_box = first_box ? temp_box : surrounding_box(temp_box, temp_box);
        first_box = false;
    }

    return temp_box;
}

double Primitive_list::pdf_value(const Point3f& o, const Vector3f& v) const {
    auto weight = 1.0 / objects.size();
    auto sum = 0.0;

    for (const auto& object : objects)
        sum += weight * object->pdf_value(o, v);

    return sum;
}

Vector3f Primitive_list::random(const Vector3f& o) const {
    auto int_size = static_cast<int>(objects.size());
    return objects[random_int(0, int_size - 1)]->random(o);
}