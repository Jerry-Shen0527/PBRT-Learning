#pragma once
#include "vec3.h"

class Medium;

class Ray {
public:
    //ray() {}
    //ray(const Point3f& origin, const Vector3f& direction, Float time = 0.0)
      //  : orig(origin), dir(direction), tm(time)
    //{}
     // Ray Public Methods
    Ray() : tMin(0.f), tMax(infinity), time(0.f) {}
    Ray(const Point3f& o, const Vector3f& d, Float time = 0.f, Float tMin_ = 0.f, Float tMax_ = infinity)
        : o(o), d(d), tMin(tMin_), tMax(tMax_), time(time) {}
    Point3f operator()(Float t) const { return o + d * t; }
    bool HasNaNs() const { return (o.HasNaNs() || d.HasNaNs() || isNaN(tMax)); }
    friend std::ostream& operator<<(std::ostream& os, const Ray& r) {
        
        os << "[o=" << r.o << ", d=" << r.d << ", tMax=" << r.tMax
            << ", time=" << r.time << "]";
        return os;
    }
    Point3f origin() const { return o; }
    Vector3f direction() const { return d; }
    Float Time() const { return time; }
    Point3f at(Float t) const {
        return o + t * d;
    }

   
    
public:
    
    Point3f o;
    Vector3f d;
    mutable Float tMin,tMax;
    Float time;
    //const Medium* medium;
    
};

using ray = Ray;
/*
class Ray {
public:
    // Ray Public Methods
    Ray() : tMax(infinity), time(0.f), medium(nullptr) {}
    Ray(const Point3f& o, const Vector3f& d, Float tMax = infinity,
        Float time = 0.f, const Medium* medium = nullptr)
        : o(o), d(d), tMax(tMax), time(time), medium(medium) {}
    Point3f operator()(Float t) const { return o + d * t; }
    bool HasNaNs() const { return (o.HasNaNs() || d.HasNaNs() || isNaN(tMax)); }
    friend std::ostream& operator<<(std::ostream& os, const Ray& r) {
        os << "[o=" << r.o << ", d=" << r.d << ", tMax=" << r.tMax
            << ", time=" << r.time << "]";
        return os;
    }

    // Ray Public Data
    Point3f o;
    Vector3f d;
    mutable Float tMax;
    Float time;
    const Medium* medium;
};
*/

class RayDifferential : public Ray {
public:
    // RayDifferential Public Methods
    RayDifferential() { hasDifferentials = false; }
    //modify---
    RayDifferential(const Point3f& o, const Vector3f& d, Float tMin = 0.f, Float tMax = Infinity,
        Float time = 0.f)
        : Ray(o, d, tMin, tMax, time) {
        hasDifferentials = false;
    }
    RayDifferential(const Ray& ray) : Ray(ray) { hasDifferentials = false; }
    bool HasNaNs() const {
        return Ray::HasNaNs() ||
            (hasDifferentials &&
                (rxOrigin.HasNaNs() || ryOrigin.HasNaNs() ||
                    rxDirection.HasNaNs() || ryDirection.HasNaNs()));
    }
    void ScaleDifferentials(Float s) {
        rxOrigin = o + (rxOrigin - o) * s;
        ryOrigin = o + (ryOrigin - o) * s;
        rxDirection = d + (rxDirection - d) * s;
        ryDirection = d + (ryDirection - d) * s;
    }
    friend std::ostream& operator<<(std::ostream& os, const RayDifferential& r) {
        os << "[ " << (Ray&)r << " has differentials: " <<
            (r.hasDifferentials ? "true" : "false") << ", xo = " << r.rxOrigin <<
            ", xd = " << r.rxDirection << ", yo = " << r.ryOrigin << ", yd = " <<
            r.ryDirection;
        return os;
    }

    // RayDifferential Public Data
    bool hasDifferentials;
    Point3f rxOrigin, ryOrigin;
    Vector3f rxDirection, ryDirection;
};


