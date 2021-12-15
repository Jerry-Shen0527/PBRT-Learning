#ifndef AABB_H
#define AABB_H

#include "rtweekend.h"
#include "vec3.h"
#include "ray.h"

//Float
class aabb {
public:
    aabb() {
        Float minNum = std::numeric_limits<Float>::lowest();
        Float maxNum = std::numeric_limits<Float>::max();
        pMin = Point3f(maxNum, maxNum, maxNum);
        pMax = Point3f(minNum, minNum, minNum);
    }
    explicit aabb(const Point3f& p) : pMin(p), pMax(p) {}
    aabb(const Point3f& p1, const Point3f& p2)
        : pMin(std::min(p1.x, p2.x), std::min(p1.y, p2.y),
            std::min(p1.z, p2.z)),
        pMax(std::max(p1.x, p2.x), std::max(p1.y, p2.y),
            std::max(p1.z, p2.z)) {}
    //aabb(const Point3f& a, const Point3f& b) { pMin = a; pMax = b; }

    Point3f min() const { return pMin; }
    Point3f max() const { return pMax; }

    bool operator==(const aabb& b) const {
        return b.pMin == pMin && b.pMax == pMax;
    }
    bool operator!=(const aabb& b) const {
        return b.pMin != pMin || b.pMax != pMax;
    }
    const Point3f& operator[](int i) const
    {
        //DCHECK(i == 0 || i == 1);
        return (i == 0) ? pMin : pMax;
    }
    Point3f& operator[](int i)
    {
        //DCHECK(i == 0 || i == 1);
        return (i == 0) ? pMin : pMax;
    }

    Point3f Corner(int corner) const {
        //DCHECK(corner >= 0 && corner < 8);
        return Point3f((*this)[(corner & 1)].x,
            (*this)[(corner & 2) ? 1 : 0].y,
            (*this)[(corner & 4) ? 1 : 0].z);
    }
    /*
    Vector3<T> Diagonal() const { return pMax - pMin; }
    T SurfaceArea() const {
        Vector3<T> d = Diagonal();
        return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    T Volume() const {
        Vector3<T> d = Diagonal();
        return d.x * d.y * d.z;
    }
    int MaximumExtent() const {
        Vector3<T> d = Diagonal();
        if (d.x > d.y && d.x > d.z)
            return 0;
        else if (d.y > d.z)
            return 1;
        else
            return 2;
    }
    Point3<T> Lerp(const Point3f& t) const {
        return Point3<T>(pbrt::Lerp(t.x, pMin.x, pMax.x),
            pbrt::Lerp(t.y, pMin.y, pMax.y),
            pbrt::Lerp(t.z, pMin.z, pMax.z));
    }
    Vector3<T> Offset(const Point3<T>& p) const {
        Vector3<T> o = p - pMin;
        if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
        if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
        if (pMax.z > pMin.z) o.z /= pMax.z - pMin.z;
        return o;
    }
    void BoundingSphere(Point3<T>* center, Float* radius) const {
        *center = (pMin + pMax) / 2;
        *radius = Inside(*center, *this) ? Distance(*center, pMax) : 0;
    }
    template <typename U>
    explicit operator Bounds3<U>() const {
        return Bounds3<U>((Point3<U>)pMin, (Point3<U>)pMax);
    }
    bool IntersectP(const Ray& ray, Float* hitt0 = nullptr,
        Float* hitt1 = nullptr) const;
    inline bool IntersectP(const Ray& ray, const Vector3f& invDir,
        const int dirIsNeg[3]) const;
    friend std::ostream& operator<<(std::ostream& os, const Bounds3<T>& b) {
        os << "[ " << b.pMin << " - " << b.pMax << " ]";
        return os;
    }*/

    bool hit(const ray& r, double t_min, double t_max) const;

    Point3f pMin, pMax;
};

inline bool aabb::hit(const ray& r, double t_min, double t_max) const {
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

inline aabb surrounding_box(aabb box0, aabb box1) {
    Point3f small(fmin(box0.min().x, box1.min().x),
        fmin(box0.min().y, box1.min().y),
        fmin(box0.min().z, box1.min().z));

    Point3f big(fmax(box0.max().x, box1.max().x),
        fmax(box0.max().y, box1.max().y),
        fmax(box0.max().z, box1.max().z));

    return aabb(small, big);
}

//Bounding unity
//template <typename T>
inline aabb Union(const aabb& b, const Point3f& p) {
    aabb ret;
    ret.pMin = Min(b.pMin, p);
    ret.pMax = Max(b.pMax, p);
    return ret;
}

//template <typename T>
inline aabb Union(const aabb& b1, const aabb& b2) {
    aabb ret;
    ret.pMin = Min(b1.pMin, b2.pMin);
    ret.pMax = Max(b1.pMax, b2.pMax);
    return ret;
}

//template <typename T>
inline aabb Intersect(const aabb& b1, const aabb& b2) {
    // Important: assign to pMin/pMax directly and don't run the Bounds2()
    // constructor, since it takes min/max of the points passed to it.  In
    // turn, that breaks returning an invalid bound for the case where we
    // intersect non-overlapping bounds (as we'd like to happen).
    aabb ret;
    ret.pMin = Max(b1.pMin, b2.pMin);
    ret.pMax = Min(b1.pMax, b2.pMax);
    return ret;
}

//template <typename T>
inline bool Overlaps(const aabb& b1, const aabb& b2) {
    bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
    bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
    bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
    return (x && y && z);
}

//template <typename T>
inline bool Inside(const Point3f& p, const aabb& b) {
    return (p.x >= b.pMin.x && p.x <= b.pMax.x && p.y >= b.pMin.y &&
        p.y <= b.pMax.y && p.z >= b.pMin.z && p.z <= b.pMax.z);
}

//template <typename T>
inline bool InsideExclusive(const Point3f& p, const aabb& b) {
    return (p.x >= b.pMin.x && p.x < b.pMax.x&& p.y >= b.pMin.y &&
        p.y < b.pMax.y&& p.z >= b.pMin.z && p.z < b.pMax.z);
}

//template <typename T, typename U>
template<typename U>
inline aabb Expand(const aabb& b, U delta) {
    return aabb(b.pMin - Vector3f(delta, delta, delta),
        b.pMax + Vector3f(delta, delta, delta));
}

// Minimum squared distance from point to box; returns zero if point is
// inside.
template <typename T>
inline Float DistanceSquared(const Point3<T>& p, const aabb& b) {
    Float dx = std::max({ Float(0), b.pMin.x - p.x, p.x - b.pMax.x });
    Float dy = std::max({ Float(0), b.pMin.y - p.y, p.y - b.pMax.y });
    Float dz = std::max({ Float(0), b.pMin.z - p.z, p.z - b.pMax.z });
    return dx * dx + dy * dy + dz * dz;
}

template <typename T>
inline Float Distance(const Point3<T>& p, const aabb& b) {
    return std::sqrt(DistanceSquared(p, b));
}

/*
template <typename T>
class Bounds3 {
public:
    // Bounds3 Public Methods
    Bounds3() {
        T minNum = std::numeric_limits<T>::lowest();
        T maxNum = std::numeric_limits<T>::max();
        pMin = Point3<T>(maxNum, maxNum, maxNum);
        pMax = Point3<T>(minNum, minNum, minNum);
    }
    explicit Bounds3(const Point3<T>& p) : pMin(p), pMax(p) {}
    Bounds3(const Point3<T>& p1, const Point3<T>& p2)
        : pMin(std::min(p1.x, p2.x), std::min(p1.y, p2.y),
            std::min(p1.z, p2.z)),
        pMax(std::max(p1.x, p2.x), std::max(p1.y, p2.y),
            std::max(p1.z, p2.z)) {}
    const Point3<T>& operator[](int i) const;
    Point3<T>& operator[](int i);
    bool operator==(const Bounds3<T>& b) const {
        return b.pMin == pMin && b.pMax == pMax;
    }
    bool operator!=(const Bounds3<T>& b) const {
        return b.pMin != pMin || b.pMax != pMax;
    }
    Point3<T> Corner(int corner) const {
        DCHECK(corner >= 0 && corner < 8);
        return Point3<T>((*this)[(corner & 1)].x,
            (*this)[(corner & 2) ? 1 : 0].y,
            (*this)[(corner & 4) ? 1 : 0].z);
    }
    Vector3<T> Diagonal() const { return pMax - pMin; }
    T SurfaceArea() const {
        Vector3<T> d = Diagonal();
        return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    T Volume() const {
        Vector3<T> d = Diagonal();
        return d.x * d.y * d.z;
    }
    int MaximumExtent() const {
        Vector3<T> d = Diagonal();
        if (d.x > d.y && d.x > d.z)
            return 0;
        else if (d.y > d.z)
            return 1;
        else
            return 2;
    }
    Point3<T> Lerp(const Point3f& t) const {
        return Point3<T>(pbrt::Lerp(t.x, pMin.x, pMax.x),
            pbrt::Lerp(t.y, pMin.y, pMax.y),
            pbrt::Lerp(t.z, pMin.z, pMax.z));
    }
    Vector3<T> Offset(const Point3<T>& p) const {
        Vector3<T> o = p - pMin;
        if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
        if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
        if (pMax.z > pMin.z) o.z /= pMax.z - pMin.z;
        return o;
    }
    void BoundingSphere(Point3<T>* center, Float* radius) const {
        *center = (pMin + pMax) / 2;
        *radius = Inside(*center, *this) ? Distance(*center, pMax) : 0;
    }
    template <typename U>
    explicit operator Bounds3<U>() const {
        return Bounds3<U>((Point3<U>)pMin, (Point3<U>)pMax);
    }
    bool IntersectP(const Ray& ray, Float* hitt0 = nullptr,
        Float* hitt1 = nullptr) const {
        Float t0 = 0, t1 = ray.tMax;
        for (int i = 0; i < 3; ++i) {
            // Update interval for _i_th bounding box slab
            Float invRayDir = 1 / ray.d[i];
            Float tNear = (pMin[i] - ray.o[i]) * invRayDir;
            Float tFar = (pMax[i] - ray.o[i]) * invRayDir;

            // Update parametric interval from slab intersection $t$ values
            if (tNear > tFar) std::swap(tNear, tFar);

            // Update _tFar_ to ensure robust ray--bounds intersection
            tFar *= 1 + 2 * gamma(3);
            t0 = tNear > t0 ? tNear : t0;
            t1 = tFar < t1 ? tFar : t1;
            if (t0 > t1) return false;
        }
        if (hitt0) *hitt0 = t0;
        if (hitt1) *hitt1 = t1;
        return true;
    }
    inline bool IntersectP(const Ray& ray, const Vector3f& invDir,
        const int dirIsNeg[3]) const {
        const Bounds3f& bounds = *this;
        // Check for ray intersection against $x$ and $y$ slabs
        Float tMin = (bounds[dirIsNeg[0]].x - ray.o.x) * invDir.x;
        Float tMax = (bounds[1 - dirIsNeg[0]].x - ray.o.x) * invDir.x;
        Float tyMin = (bounds[dirIsNeg[1]].y - ray.o.y) * invDir.y;
        Float tyMax = (bounds[1 - dirIsNeg[1]].y - ray.o.y) * invDir.y;

        // Update _tMax_ and _tyMax_ to ensure robust bounds intersection
        tMax *= 1 + 2 * gamma(3);
        tyMax *= 1 + 2 * gamma(3);
        if (tMin > tyMax || tyMin > tMax) return false;
        if (tyMin > tMin) tMin = tyMin;
        if (tyMax < tMax) tMax = tyMax;

        // Check for ray intersection against $z$ slab
        Float tzMin = (bounds[dirIsNeg[2]].z - ray.o.z) * invDir.z;
        Float tzMax = (bounds[1 - dirIsNeg[2]].z - ray.o.z) * invDir.z;

        // Update _tzMax_ to ensure robust bounds intersection
        tzMax *= 1 + 2 * gamma(3);
        if (tMin > tzMax || tzMin > tMax) return false;
        if (tzMin > tMin) tMin = tzMin;
        if (tzMax < tMax) tMax = tzMax;
        return (tMin < ray.tMax) && (tMax > 0);
    }
    friend std::ostream& operator<<(std::ostream& os, const Bounds3<T>& b) {
        os << "[ " << b.pMin << " - " << b.pMax << " ]";
        return os;
    }

    // Bounds3 Public Data
    Point3<T> pMin, pMax;
};

template <typename T>
Bounds3<T> Union(const Bounds3<T> &b, const Point3<T> &p) {
    Bounds3<T> ret;
    ret.pMin = Min(b.pMin, p);
    ret.pMax = Max(b.pMax, p);
    return ret;
}

template <typename T>
Bounds3<T> Union(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    Bounds3<T> ret;
    ret.pMin = Min(b1.pMin, b2.pMin);
    ret.pMax = Max(b1.pMax, b2.pMax);
    return ret;
}

template <typename T>
Bounds3<T> Intersect(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    // Important: assign to pMin/pMax directly and don't run the Bounds2()
    // constructor, since it takes min/max of the points passed to it.  In
    // turn, that breaks returning an invalid bound for the case where we
    // intersect non-overlapping bounds (as we'd like to happen).
    Bounds3<T> ret;
    ret.pMin = Max(b1.pMin, b2.pMin);
    ret.pMax = Min(b1.pMax, b2.pMax);
    return ret;
}

template <typename T>
bool Overlaps(const Bounds3<T> &b1, const Bounds3<T> &b2) {
    bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
    bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
    bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
    return (x && y && z);
}

template <typename T>
bool Inside(const Point3<T> &p, const Bounds3<T> &b) {
    return (p.x >= b.pMin.x && p.x <= b.pMax.x && p.y >= b.pMin.y &&
            p.y <= b.pMax.y && p.z >= b.pMin.z && p.z <= b.pMax.z);
}

template <typename T>
bool InsideExclusive(const Point3<T> &p, const Bounds3<T> &b) {
    return (p.x >= b.pMin.x && p.x < b.pMax.x && p.y >= b.pMin.y &&
            p.y < b.pMax.y && p.z >= b.pMin.z && p.z < b.pMax.z);
}

template <typename T, typename U>
inline Bounds3<T> Expand(const Bounds3<T> &b, U delta) {
    return Bounds3<T>(b.pMin - Vector3<T>(delta, delta, delta),
                      b.pMax + Vector3<T>(delta, delta, delta));
}

// Minimum squared distance from point to box; returns zero if point is
// inside.
template <typename T, typename U>
inline Float DistanceSquared(const Point3<T> &p, const Bounds3<U> &b) {
    Float dx = std::max({Float(0), b.pMin.x - p.x, p.x - b.pMax.x});
    Float dy = std::max({Float(0), b.pMin.y - p.y, p.y - b.pMax.y});
    Float dz = std::max({Float(0), b.pMin.z - p.z, p.z - b.pMax.z});
    return dx * dx + dy * dy + dz * dz;
}

template <typename T, typename U>
inline Float Distance(const Point3<T> &p, const Bounds3<U> &b) {
    return std::sqrt(DistanceSquared(p, b));
}
*/
#endif