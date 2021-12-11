#ifndef AABB_H
#define AABB_H

#include"ray.h"


class aabb {
public:
    aabb() {
        pMin = point3(-Infinity, -Infinity, -Infinity);
        pMax = point3(Infinity, Infinity, Infinity);
    }
    aabb(const point3& p):pMin(p),pMax(p){}
    aabb(const point3& p1, const point3& p2)  :pMin(std::min(p1.x, p2.x), std::min(p1.y, p2.y),std::min(p1.z, p2.z)),
        pMax(std::max(p1.x, p2.x), std::max(p1.y, p2.y),std::max(p1.z, p2.z)) { }

    inline point3 const operator[](int i) const{
        IN_RANGE((i == 0 || i == 1));
        return (i == 0) ? pMin : pMax;
    }
    inline point3 &operator[](int i) {
        IN_RANGE((i == 0 || i == 1));
        return (i == 0) ? pMin : pMax;
    }

    point3 min() const { return pMin; }
    point3 max() const { return pMax; }

    inline bool hit(const ray& r, Float t_min, Float t_max) const {
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
    bool Overlaps(const aabb& b) const
    {
        bool x = (pMax.x >= b.pMin.x) && (pMin.x <= b.pMax.x);
        bool y = (pMax.y >= b.pMin.y) && (pMin.y <= b.pMax.y);
        bool z = (pMax.z >= b.pMin.z) && (pMin.z <= b.pMax.z);
        return x && y && z;
    }
    bool Inside(const point3& b) const
    {
        bool x = (pMax.x >= b.x) && (pMin.x <= b.x);
        bool y = (pMax.y >= b.y) && (pMin.y <= b.y);
        bool z = (pMax.z >= b.z) && (pMin.z <= b.z);
        return x && y && z;
    }
    void Expand(Float delta)
    {
        pMin -= vec3(delta, delta, delta);
        pMax += vec3(delta, delta, delta);
    }
    Float SurfaceArea() const
    {
        vec3 d = pMax - pMin;
        return 2.f * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    Float Volume()const
    {
        vec3 d = pMax - pMin;
        return d.x * d.y * d.z;
    }

    /*Return longest axis*/
    int MaximumExtent() const
    {
        vec3 diag = pMax - pMin;
        if (diag.x > diag.y && diag.x > diag.z)
            return 0;
        else if (diag.y > diag.z)
            return 1;
        else
            return 2;
    }

    Point3f Corner(int corner) const {
        IN_RANGE((corner >= 0 && corner < 8));
        return Point3f((*this)[(corner & 1)].x,
            (*this)[(corner & 2) ? 1 : 0].y,
            (*this)[(corner & 4) ? 1 : 0].z);
    }

    point3 Lerp(Float tx, Float ty, Float tz) const
    {
        return point3(::Lerp(tx, pMin.x, pMax.y), ::Lerp(ty, pMin.y, pMax.y), ::Lerp(tz, pMin.z, pMax.z));
    }
    /*Return relative position in box*/
    vec3 Offset(const point3& p) const
    {
        return vec3((p.x - pMin.x) / (pMax.x - pMin.x),
                    (p.y - pMin.y) / (pMax.y - pMin.y),
                    (p.z - pMin.z) / (pMax.z - pMin.z));
    }
    void BoundingSphere(point3& c, float& rad) const
    {
        c = .5f * pMin + .5f * pMax;
        rad=Inside(c)?Distance(c,pMax):0.f;
    }
    /*计算是否在ray的(0,tMax)范围*/
    bool aabb::IntersectP(const Ray& ray, Float* hitt0, Float* hitt1) const;
    /*通过预计算而不反复交换，但不计算是否符合ray范围*/
    bool aabb::IntersectP(const Ray& ray, const Vector3f& invDir, const int dirIsNeg[3]) const;


    point3 pMin;
    point3 pMax;
};

inline aabb Union(aabb box0, aabb box1) {
    point3 small(fmin(box0.min().x, box1.min().x),
        fmin(box0.min().y, box1.min().y),
        fmin(box0.min().z, box1.min().z));

    point3 big(fmax(box0.max().x, box1.max().x),
        fmax(box0.max().y, box1.max().y),
        fmax(box0.max().z, box1.max().z));

    return aabb(small, big);
}


using  Bounds3f = aabb;
#endif

