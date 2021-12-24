#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_SHAPES_SHAPE_LIST_H
#define PBRT_SHAPES_SHAPE_LIST_H

// shapes/triangle.h*
#include "shape.h"
//#include "stats.h"
#include <map>
#include<string>

namespace pbrt {

    class Shape_list :public Shape {
    public:
        //Shape_list() {}
        Shape_list(const std::vector<shared_ptr<Shape>>& Objects,
            shared_ptr<Transform> ObjectToWorld , shared_ptr<Transform> WorldToObject ,
            bool reverseOrientation=false)
            : Shape(ObjectToWorld, WorldToObject, reverseOrientation) {
            objects = Objects;
        }

        void clear() { objects.clear(); }
        void add(shared_ptr<Shape> object) { objects.push_back(object); }
        Bounds3f ObjectBound() const;
        Bounds3f WorldBound() const;
        bool Intersect(const Ray& ray, Float* tHit, SurfaceInteraction* isect,
            bool testAlphaTexture = true) const;
        bool IntersectP(const Ray& ray, bool testAlphaTexture = true) const;
        Float Area() const;

    public:
        std::vector<shared_ptr<Shape>> objects;
    };
    

}  // namespace pbrt

#endif  // PBRT_SHAPES_SHAPE_LIST_H