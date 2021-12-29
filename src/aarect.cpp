#include "aarect.h"

bool xy_rect::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    auto t = (k - r.origin().z) / r.direction().z;
    if (t < t_min || t > t_max)
        return false;
    auto x = r.origin().x + t * r.direction().x;
    auto y = r.origin().y + t * r.direction().y;
    if (x < x0 || x > x1 || y < y0 || y > y1)
        return false;
    rec.u = (x - x0) / (x1 - x0);
    rec.v = (y - y0) / (y1 - y0);
    rec.t = t;
    auto outward_normal = Normal3f(0, 0, 1);
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mp;
    rec.p = r.at(t);
    return true;
}

bool xz_rect::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    auto t = (k - r.origin().y) / r.direction().y;
    if (t < t_min || t > t_max)
        return false;
    auto x = r.origin().x + t * r.direction().x;
    auto z = r.origin().z + t * r.direction().z;
    if (x < x0 || x > x1 || z < z0 || z > z1)
        return false;
    rec.u = (x - x0) / (x1 - x0);
    rec.v = (z - z0) / (z1 - z0);
    rec.t = t;
    auto outward_normal = Normal3f(0, 1, 0);
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mp;
    rec.p = r.at(t);
    return true;
}

bool yz_rect::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    auto t = (k - r.origin().x) / r.direction().x;
    if (t < t_min || t > t_max)
        return false;
    auto y = r.origin().y + t * r.direction().y;
    auto z = r.origin().z + t * r.direction().z;
    if (y < y0 || y > y1 || z < z0 || z > z1)
        return false;
    rec.u = (y - y0) / (y1 - y0);
    rec.v = (z - z0) / (z1 - z0);
    rec.t = t;
    auto outward_normal = Normal3f(1, 0, 0);
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mp;
    rec.p = r.at(t);
    return true;
}

namespace pbrt {
    Bounds3f XY_Rect::ObjectBound() const{
        return Bounds3f(Point3f(x0, y0, z),
            Point3f(x1, y1, z));
    }

    bool XY_Rect::Intersect(const Ray& r, Float* tHit, SurfaceInteraction* isect,
        bool testAlphaTexture) const {
        //ProfilePhase p(Prof::ShapeIntersect);
        // Transform _Ray_ to object space
        Vector3f oErr, dErr;
        Ray ray = (*WorldToObject)(r, &oErr, &dErr);

        // Compute plane intersection for XY_Rect

        // Reject XY_Rect intersections for rays parallel to the XY_Rect's plane
        if (ray.d.z == 0) return false;
        Float tShapeHit = (z - ray.o.z) / ray.d.z;
        if (tShapeHit <= 0 || tShapeHit >= ray.tMax) return false;

        // See if hit point is inside disk radii and $\phimax$
        Point3f pHit = ray(tShapeHit);
        if (pHit.x < x0 || pHit.x > x1 || pHit.y < y0 || pHit.y > y1)
            return false;

        // Find parametric representation of XY_Rect hit
        Float u = (pHit.x - x0) / (x1 - x0);
        Float v = (pHit.y - y0) / (y1 - y0);
        Vector3f dpdu(x1 - x0, 0, 0);
        Vector3f dpdv(0, y1 - y0, 0);
        Normal3f dndu(0, 0, 0), dndv(0, 0, 0);

        // Refine disk intersection point
        pHit.z = z;

        // Compute error bounds for disk intersection
        Vector3f pError(0, 0, 0);

        // Initialize _SurfaceInteraction_ from parametric information
        *isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(u, v),
            -ray.d, dpdu, dpdv, dndu, dndv,
            ray.time, this));

        // Update _tHit_ for quadric intersection
        *tHit = (Float)tShapeHit;
        return true;

    }

    bool XY_Rect::IntersectP(const Ray& r, bool testAlphaTexture) const {
        //ProfilePhase p(Prof::ShapeIntersectP);
        // Transform _Ray_ to object space
        Vector3f oErr, dErr;
        Ray ray = (*WorldToObject)(r, &oErr, &dErr);

        // Compute plane intersection for XY_Rect

        // Reject XY_Rect intersections for rays parallel to the XY_Rect's plane
        if (ray.d.z == 0) return false;
        Float tShapeHit = (z - ray.o.z) / ray.d.z;
        if (tShapeHit <= 0 || tShapeHit >= ray.tMax) return false;

        // See if hit point is inside disk radii and $\phimax$
        Point3f pHit = ray(tShapeHit);
        if (pHit.x < x0 || pHit.x > x1 || pHit.y < y0 || pHit.y > y1)
            return false;

        return true;
    }

    Float XY_Rect::Area() const {
        return (x1 - x0) * (y1 - y0);
    }

    //Interaction Disk::Sample(const Point2f& u, Float* pdf) const {
    //    Point2f pd = ConcentricSampleDisk(u);
    //    Point3f pObj(pd.x * radius, pd.y * radius, height);
    //    Interaction it;
    //    it.n = Normalize((*ObjectToWorld)(Normal3f(0, 0, 1)));
    //    if (reverseOrientation) it.n *= -1;
    //    it.p = (*ObjectToWorld)(pObj, Vector3f(0, 0, 0), &it.pError);
    //    *pdf = 1 / Area();
    //    return it;
    //}


    Bounds3f YZ_Rect::ObjectBound() const {
        return Bounds3f(Point3f(y0, z0, x),
            Point3f(y1, z1, x));
    }

    bool YZ_Rect::Intersect(const Ray& r, Float* tHit, SurfaceInteraction* isect,
        bool testAlphaTexture) const {
        //ProfilePhase p(Prof::ShapeIntersect);
        // Transform _Ray_ to object space
        Vector3f oErr, dErr;
        Ray ray = (*WorldToObject)(r, &oErr, &dErr);

        // Compute plane intersection for XY_Rect

        // Reject XY_Rect intersections for rays parallel to the XY_Rect's plane
        if (ray.d.x == 0) return false;
        Float tShapeHit = (x - ray.o.x) / ray.d.x;
        if (tShapeHit <= 0 || tShapeHit >= ray.tMax) return false;

        // See if hit point is inside disk radii and $\phimax$
        Point3f pHit = ray(tShapeHit);
        if (pHit.y < y0 || pHit.y > y1 || pHit.z < z0 || pHit.z > z1)
            return false;

        // Find parametric representation of XY_Rect hit
        Float u = (pHit.y - y0) / (y1 - y0);
        Float v = (pHit.z - z0) / (z1 - z0);
        Vector3f dpdu(0, y1 - y0, 0);
        Vector3f dpdv(0, 0, z1 - z0);
        Normal3f dndu(0, 0, 0), dndv(0, 0, 0);

        // Refine disk intersection point
        pHit.x = x;

        // Compute error bounds for disk intersection
        Vector3f pError(0, 0, 0);

        // Initialize _SurfaceInteraction_ from parametric information
        *isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(u, v),
            -ray.d, dpdu, dpdv, dndu, dndv,
            ray.time, this));

        // Update _tHit_ for quadric intersection
        *tHit = (Float)tShapeHit;
        return true;

    }

    bool YZ_Rect::IntersectP(const Ray& r, bool testAlphaTexture) const {
        //ProfilePhase p(Prof::ShapeIntersectP);
        // Transform _Ray_ to object space
        Vector3f oErr, dErr;
        Ray ray = (*WorldToObject)(r, &oErr, &dErr);

        // Compute plane intersection for XY_Rect

        // Reject XY_Rect intersections for rays parallel to the XY_Rect's plane
        if (ray.d.x == 0) return false;
        Float tShapeHit = (x - ray.o.x) / ray.d.x;
        if (tShapeHit <= 0 || tShapeHit >= ray.tMax) return false;

        // See if hit point is inside disk radii and $\phimax$
        Point3f pHit = ray(tShapeHit);
        if (pHit.y < y0 || pHit.y > y1 || pHit.z < z0 || pHit.z > z1)
            return false;

        return true;
    }

    Float YZ_Rect::Area() const {
        return (y1 - y0) * (z1 - z0);
    }

    //ZX
    Bounds3f ZX_Rect::ObjectBound() const {
        return Bounds3f(Point3f(z0, x0, y),
            Point3f(z1, x1, y));
    }

    bool ZX_Rect::Intersect(const Ray& r, Float* tHit, SurfaceInteraction* isect,
        bool testAlphaTexture) const {
        //ProfilePhase p(Prof::ShapeIntersect);
        // Transform _Ray_ to object space
        Vector3f oErr, dErr;
        Ray ray = (*WorldToObject)(r, &oErr, &dErr);

        // Compute plane intersection for XY_Rect

        // Reject XY_Rect intersections for rays parallel to the XY_Rect's plane
        if (ray.d.y == 0) return false;
        Float tShapeHit = (y - ray.o.y) / ray.d.y;
        if (tShapeHit <= 0 || tShapeHit >= ray.tMax) return false;

        // See if hit point is inside disk radii and $\phimax$
        Point3f pHit = ray(tShapeHit);
        if (pHit.z < z0 || pHit.z > z1 || pHit.x < x0 || pHit.x > x1)
            return false;

        // Find parametric representation of XY_Rect hit
        Float u = (pHit.z - z0) / (z1 - z0);
        Float v = (pHit.x - x0) / (x1 - x0);
        Vector3f dpdu(0, 0, z1 - z0);
        Vector3f dpdv(x1 - x0, 0, 0);
        Normal3f dndu(0, 0, 0), dndv(0, 0, 0);

        // Refine disk intersection point
        pHit.y = y;

        // Compute error bounds for disk intersection
        Vector3f pError(0, 0, 0);

        // Initialize _SurfaceInteraction_ from parametric information
        *isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(u, v),
            -ray.d, dpdu, dpdv, dndu, dndv,
            ray.time, this));

        // Update _tHit_ for quadric intersection
        *tHit = (Float)tShapeHit;
        return true;

    }

    bool ZX_Rect::IntersectP(const Ray& r, bool testAlphaTexture) const {
        //ProfilePhase p(Prof::ShapeIntersectP);
        // Transform _Ray_ to object space
        Vector3f oErr, dErr;
        Ray ray = (*WorldToObject)(r, &oErr, &dErr);

        // Compute plane intersection for XY_Rect

        // Reject XY_Rect intersections for rays parallel to the XY_Rect's plane
        if (ray.d.y == 0) return false;
        Float tShapeHit = (y - ray.o.y) / ray.d.y;
        if (tShapeHit <= 0 || tShapeHit >= ray.tMax) return false;

        // See if hit point is inside disk radii and $\phimax$
        Point3f pHit = ray(tShapeHit);
        if (pHit.z < z0 || pHit.z > z1 || pHit.x < x0 || pHit.x > x1)
            return false;

        return true;
    }

    Float ZX_Rect::Area() const {
        return (z1 - z0) * (x1 - x0);
    }


}