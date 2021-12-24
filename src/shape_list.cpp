#include "shape_list.h"

namespace pbrt {
    Bounds3f Shape_list::ObjectBound() const{
        Bounds3f objBound;
        for (auto& shape : objects) {
            objBound = Union(objBound, shape->ObjectBound());
        }
        return objBound;
    }
    Bounds3f Shape_list::WorldBound() const{
        Bounds3f worBound;
        for (auto& shape : objects) {
            worBound = Union(worBound, shape->WorldBound());
        }
        return worBound;
    }
    bool Shape_list::Intersect(const Ray& ray, Float* tHit, SurfaceInteraction* isect,
        bool testAlphaTexture) const{
        bool hit_list = false;
        for (auto& shape : objects) {
            if (shape->IntersectP(ray, testAlphaTexture))
            {
                hit_list = true;
                shape->Intersect(ray, tHit, isect, testAlphaTexture);
                ray.tMax = *tHit;
            }
        }
        if (hit_list)
            return true;
        return false;
    }

    bool Shape_list::IntersectP(const Ray& ray, bool testAlphaTexture) const{
        for (auto& shape : objects) {
            if (shape->IntersectP(ray, testAlphaTexture))
            return true;
        }
        return false;
    }
    Float Shape_list::Area() const{
        Float are = 0.;
        for (auto& shape : objects) {
            are += shape->Area();
        }
        return are;
    }
}