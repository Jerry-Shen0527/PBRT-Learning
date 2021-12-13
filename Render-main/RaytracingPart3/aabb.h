#ifndef AABB_H
#define AABB_H

#include "rtweekend.h"
#include "efloat.h"
class aabb {
public:
    aabb() {}
    aabb(const Point3f& a, const Point3f& b) { minimum = a; maximum = b; }

    bool operator==(const aabb& b) const {
        return b.minimum == minimum && b.maximum == maximum;
    }
    bool operator!=(const aabb& b) const {
        return b.minimum != minimum || b.maximum != maximum;
    }
    
    Vector3f Diagonal() const { return maximum - minimum; }
    double SurfaceArea() const {
        Vector3f d = Diagonal();
        return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    double Volume() const {
        Vector3f d = Diagonal();
        return d.x * d.y * d.z;
    }
    int MaximumExtent() const {
        Vector3f d = Diagonal();
        if (d.x > d.y && d.x > d.z)
            return 0;
        else if (d.y > d.z)
            return 1;
        else
            return 2;
    }

    Vector3f Offset(const Point3f& p) const {
        Vector3f o = p - minimum;
        if (maximum.x > minimum.x) o.x /= maximum.x - minimum.x;
        if (maximum.y > minimum.y) o.y /= maximum.y - minimum.y;
        if (maximum.z > minimum.z) o.z /= maximum.z - minimum.z;
        return o;
    }

    Point3f min() const { return minimum; }
    Point3f max() const { return maximum; }

    bool hit(const ray& r, double t_min, double t_max) const {
        for (int a = 0; a < 3; a++) {
            auto t0 = fmin((minimum[a] - r.origin()[a]) / r.direction()[a],(maximum[a] - r.origin()[a]) / r.direction()[a]);
            auto t1 = fmax((minimum[a] - r.origin()[a]) / r.direction()[a],(maximum[a] - r.origin()[a]) / r.direction()[a]);
            t_min = fmax(t0, t_min);
            t_max = fmin(t1, t_max);
            if (t_max <= t_min)
                return false;
        }
        return true;
    }

    inline const Point3f& operator[](int i) const {
        assert(i < 2);
        return (i == 0) ? minimum : maximum;
    }

    template <typename T>
    inline Point3f& operator[](int i) {
        assert(i < 2);
        return (i == 0) ? minimum : maximum;
    }

    inline bool IntersectP(const ray& ray, double t_min, double t_max,const Vector3f& invDir,
        const int dirIsNeg[3]) const {
        const aabb& bounds = *this;
        // Check for ray intersection against $x$ and $y$ slabs
        Float tMin = (bounds[dirIsNeg[0]].x - ray.origin().x) * invDir.x;
        Float tMax = (bounds[1 - dirIsNeg[0]].x - ray.origin().x) * invDir.x;
        Float tyMin = (bounds[dirIsNeg[1]].y - ray.origin().y) * invDir.y;
        Float tyMax = (bounds[1 - dirIsNeg[1]].y - ray.origin().y) * invDir.y;

        // Update _tMax_ and _tyMax_ to ensure robust bounds intersection
        tMax *= 1 + 2 * gamma(3);
        tyMax *= 1 + 2 * gamma(3);
        if (tMin > tyMax || tyMin > tMax) return false;
        if (tyMin > tMin) tMin = tyMin;
        if (tyMax < tMax) tMax = tyMax;

        // Check for ray intersection against $z$ slab
        Float tzMin = (bounds[dirIsNeg[2]].z - ray.origin().z) * invDir.z;
        Float tzMax = (bounds[1 - dirIsNeg[2]].z - ray.origin().z) * invDir.z;

        // Update _tzMax_ to ensure robust bounds intersection
        tzMax *= 1 + 2 * gamma(3);
        if (tMin > tzMax || tzMin > tMax) return false;
        if (tzMin > tMin) tMin = tzMin;
        if (tzMax < tMax) tMax = tzMax;
        return (tMin < t_max) && (tMax > 0);
    }

    inline bool IntersectP(const ray& ray, double t_min, double t_max,Float* hitt0,
        Float* hitt1) const {
        Float t0 = 0, t1 = t_max;
        for (int i = 0; i < 3; ++i) {
            // Update interval for _i_th bounding box slab
            Float invRayDir = 1 / ray.direction()[i];
            Float tNear = (minimum[i] - ray.origin()[i]) * invRayDir;
            Float tFar = (maximum[i] - ray.origin()[i]) * invRayDir;

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

    Point3f minimum;
    Point3f maximum;
};

aabb surrounding_box(aabb box0, aabb box1);

inline aabb Union(const aabb& b, const Point3f& p) {
    aabb ret;
    ret.minimum = Min(b.minimum, p);
    ret.maximum = Max(b.maximum, p);
    return ret;
}

inline aabb Union(const aabb& b1, const aabb& b2) {
    aabb ret;
    ret.minimum = Min(b1.minimum, b2.minimum);
    ret.maximum = Max(b1.maximum, b2.maximum);
    return ret;
}

#endif
