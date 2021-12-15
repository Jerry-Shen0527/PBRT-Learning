#pragma once
#include "vec3.h"

class Medium;

class ray {
public:
    ray() {}
    ray(const Point3f& origin, const Vector3f& direction, Float time = 0.0)
        : orig(origin), dir(direction), tm(time)
    {}

    Point3f origin() const { return orig; }
    Vector3f direction() const { return dir; }
    double time() const { return tm; }

    Point3f operator()(Float t) const { return orig + t * dir; }
    Point3f at(Float t) const {
        return orig + t * dir;
    }

    bool HasNaNs() const { return (orig.HasNaNs() || dir.HasNaNs() ); }
    //bool HasNaNs() const { return (orig.HasNaNs() || dir.HasNaNs() || isNaN(tMax)); }
    friend std::ostream& operator<<(std::ostream& os, const ray& r) {
        os << "[orig=" << r.orig << ", dir=" << r.dir << ", time=" << r.tm << "]";
        return os;
    }
public:
    Point3f orig;
    Vector3f dir;
    //mutable Float tMax;
    Float tm;//Ray with time information
};


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

class RayDifferential : public ray {
public:
    // RayDifferential Public Methods
    RayDifferential() { hasDifferentials = false; }
    RayDifferential(const Point3f& o, const Vector3f& d, Float tMax = infinity,
        Float time = 0.f, const Medium* medium = nullptr)
        : ray(o, d,  time) {
        hasDifferentials = false;
    }
    RayDifferential(const ray& ray) : ray(ray) { hasDifferentials = false; }
    bool HasNaNs() const {
        return ray::HasNaNs() ||
            (hasDifferentials &&
                (rxOrigin.HasNaNs() || ryOrigin.HasNaNs() ||
                    rxDirection.HasNaNs() || ryDirection.HasNaNs()));
    }
    void ScaleDifferentials(Float s) {
        rxOrigin = orig + (rxOrigin - orig) * s;
        ryOrigin = orig + (ryOrigin - orig) * s;
        rxDirection = dir + (rxDirection - dir) * s;
        ryDirection = dir + (ryDirection - dir) * s;
    }
    friend std::ostream& operator<<(std::ostream& os, const RayDifferential& r) {
        os << "[ " << (ray&)r << " has differentials: " <<
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

