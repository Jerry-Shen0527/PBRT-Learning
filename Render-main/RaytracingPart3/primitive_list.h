#pragma once
#include "primitive.h"

#include <memory>
#include <vector>

using std::shared_ptr;
using std::make_shared;
class Primitive_list : public Primitive
{
public:
    Primitive_list() {}
    Primitive_list(shared_ptr<Primitive> object) { add(object); }

    void clear() { objects.clear(); }
    void add(shared_ptr<Primitive> object) { objects.push_back(object); }
    void add(vector<shared_ptr<Primitive>> object) { 
        for(int i=0;i<object.size();i++)
        objects.push_back(object[i]); }

    virtual bool Intersect(const ray& r, double t_min, double t_max, SurfaceInteraction& inter_record) const override;
    virtual bool IntersectP(const ray& r, double t_min, double t_max) const override;

    virtual aabb bounding_box(
        ) const override;

    double pdf_value(const Point3f& o, const Vector3f& v) const;
    Vector3f random(const Vector3f& o) const;

public:
    std::vector<shared_ptr<Primitive>> objects;
};

