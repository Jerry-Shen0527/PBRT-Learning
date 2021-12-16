#include "sphere.h"
#include "onb.h"

bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    Vector3f oc = r.origin() - center;
    auto a = r.direction().LengthSquared();
    auto half_b = Dot(oc, r.direction());
    auto c = oc.LengthSquared() - radius * radius;

    auto discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return false;
    auto sqrtd = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    auto root = (-half_b - sqrtd) / a;
    if (root < t_min || t_max < root) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || t_max < root)
            return false;
    }

    rec.t = root;
    rec.p = r.at(rec.t);
    Normal3f outward_normal = Normal3f((rec.p - center) / radius);
    rec.set_face_normal(r, outward_normal);
    get_sphere_uv(Vector3f(outward_normal), rec.u, rec.v);
    rec.mat_ptr = mat_ptr;

    return true;
}

bool sphere::bounding_box(double time0, double time1, aabb& output_box) const {
    output_box = aabb(
        center - Vector3f(radius, radius, radius),
        center + Vector3f(radius, radius, radius));
    return true;
}

void sphere::get_sphere_uv(const Vector3f& p, double& u, double& v) {
    // p: a given point on the sphere of radius one, centered at the origin.
    // u: returned value [0,1] of angle around the Y axis from X=-1.
    // v: returned value [0,1] of angle from Y=-1 to Y=+1.
    //     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
    //     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
    //     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>

    auto theta = acos(-p.y);
    auto phi = atan2(-p.z, p.x) + pi;

    u = phi / (2 * pi);
    v = theta / pi;
}

double sphere::pdf_value(const Point3f& o, const Vector3f& v) const {
    hit_record rec;
    if (!this->hit(ray(o, v), 0.001, infinity, rec))
        return 0;

    auto cos_theta_max = sqrt(1 - radius * radius / (center - o).LengthSquared());
    auto solid_angle = 2 * pi * (1 - cos_theta_max);

    return  1 / solid_angle;
}



Vector3f sphere::random(const Point3f& o) const {
    Vector3f direction = center - o;
    auto distance_squared = direction.LengthSquared();
    onb uvw;
    uvw.build_from_w(direction);
    return uvw.local(random_to_sphere(radius, distance_squared));
}


//#include "shapes/sphere.h"
//#include "sampling.h"
//#include "paramset.h"
#include "efloat.h"
//#include "stats.h"

namespace pbrt {

    // Sphere Method Definitions
    Bounds3f Sphere::ObjectBound() const {
        return Bounds3f(Point3f(-radius, -radius, zMin),
            Point3f(radius, radius, zMax));
    }

    bool Sphere::Intersect(const Ray& r, Float* tHit, SurfaceInteraction* isect,
        bool testAlphaTexture) const {
        //ProfilePhase p(Prof::ShapeIntersect);
        Float phi;
        Point3f pHit;
        // Transform _Ray_ to object space
        Vector3f oErr, dErr;
        Ray ray = (*WorldToObject)(r, &oErr, &dErr);

        // Compute quadratic sphere coefficients

        // Initialize _EFloat_ ray coordinate values
        EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
        EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
        EFloat a = dx * dx + dy * dy + dz * dz;
        EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
        EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

        // Solve quadratic equation for _t_ values
        EFloat t0, t1;
        if (!Quadratic(a, b, c, &t0, &t1)) return false;

        // Check quadric shape _t0_ and _t1_ for nearest intersection
        //if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) return false;
        if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= ray.tMin) return false;
        EFloat tShapeHit = t0;
        if (tShapeHit.LowerBound() <= 0) {
            tShapeHit = t1;
            if (tShapeHit.UpperBound() > ray.tMax) return false;
        }

        // Compute sphere hit position and $\phi$
        pHit = ray((Float)tShapeHit);

        // Refine sphere intersection point
        pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
        if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) phi += 2 * Pi;

        // Test sphere intersection against clipping parameters
        if ((zMin > -radius && pHit.z < zMin) || (zMax < radius && pHit.z > zMax) ||
            phi > phiMax) {
            if (tShapeHit == t1) return false;
            if (t1.UpperBound() > ray.tMax) return false;
            tShapeHit = t1;
            // Compute sphere hit position and $\phi$
            pHit = ray((Float)tShapeHit);

            // Refine sphere intersection point
            pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
            if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
            phi = std::atan2(pHit.y, pHit.x);
            if (phi < 0) phi += 2 * Pi;
            if ((zMin > -radius && pHit.z < zMin) ||
                (zMax < radius && pHit.z > zMax) || phi > phiMax)
                return false;
        }

        // Find parametric representation of sphere hit
        Float u = phi / phiMax;
        Float theta = std::acos(Clamp(pHit.z / radius, -1, 1));
        Float v = (theta - thetaMin) / (thetaMax - thetaMin);

        // Compute sphere $\dpdu$ and $\dpdv$
        Float zRadius = std::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
        Float invZRadius = 1 / zRadius;
        Float cosPhi = pHit.x * invZRadius;
        Float sinPhi = pHit.y * invZRadius;
        Vector3f dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
        Vector3f dpdv =
            (thetaMax - thetaMin) *
            Vector3f(pHit.z * cosPhi, pHit.z * sinPhi, -radius * std::sin(theta));

        // Compute sphere $\dndu$ and $\dndv$
        Vector3f d2Pduu = -phiMax * phiMax * Vector3f(pHit.x, pHit.y, 0);
        Vector3f d2Pduv =
            (thetaMax - thetaMin) * pHit.z * phiMax * Vector3f(-sinPhi, cosPhi, 0.);
        Vector3f d2Pdvv = -(thetaMax - thetaMin) * (thetaMax - thetaMin) *
            Vector3f(pHit.x, pHit.y, pHit.z);

        // Compute coefficients for fundamental forms
        Float E = Dot(dpdu, dpdu);
        Float F = Dot(dpdu, dpdv);
        Float G = Dot(dpdv, dpdv);
        Vector3f N = Normalize(Cross(dpdu, dpdv));
        Float e = Dot(N, d2Pduu);
        Float f = Dot(N, d2Pduv);
        Float g = Dot(N, d2Pdvv);

        // Compute $\dndu$ and $\dndv$ from fundamental form coefficients
        Float invEGF2 = 1 / (E * G - F * F);
        Normal3f dndu = Normal3f((f * F - e * G) * invEGF2 * dpdu +
            (e * F - f * E) * invEGF2 * dpdv);
        Normal3f dndv = Normal3f((g * F - f * G) * invEGF2 * dpdu +
            (f * F - g * E) * invEGF2 * dpdv);

        // Compute error bounds for sphere intersection
        Vector3f pError = gamma(5) * Abs((Vector3f)pHit);

        // Initialize _SurfaceInteraction_ from parametric information
        *isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(u, v),
            -ray.d, dpdu, dpdv, dndu, dndv,
            ray.time, this));

        // Update _tHit_ for quadric intersection
        *tHit = (Float)tShapeHit;
        return true;
    }

    bool Sphere::IntersectP(const Ray& r, bool testAlphaTexture) const {
        //ProfilePhase p(Prof::ShapeIntersectP);
        Float phi;
        Point3f pHit;
        // Transform _Ray_ to object space
        Vector3f oErr, dErr;
        Ray ray = (*WorldToObject)(r, &oErr, &dErr);

        // Compute quadratic sphere coefficients

        // Initialize _EFloat_ ray coordinate values
        EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
        EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
        EFloat a = dx * dx + dy * dy + dz * dz;
        EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
        EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

        // Solve quadratic equation for _t_ values
        EFloat t0, t1;
        if (!Quadratic(a, b, c, &t0, &t1)) return false;

        // Check quadric shape _t0_ and _t1_ for nearest intersection
        if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) return false;
        EFloat tShapeHit = t0;
        if (tShapeHit.LowerBound() <= 0) {
            tShapeHit = t1;
            if (tShapeHit.UpperBound() > ray.tMax) return false;
        }

        // Compute sphere hit position and $\phi$
        pHit = ray((Float)tShapeHit);

        // Refine sphere intersection point
        pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
        if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) phi += 2 * Pi;

        // Test sphere intersection against clipping parameters
        if ((zMin > -radius && pHit.z < zMin) || (zMax < radius && pHit.z > zMax) ||
            phi > phiMax) {
            if (tShapeHit == t1) return false;
            if (t1.UpperBound() > ray.tMax) return false;
            tShapeHit = t1;
            // Compute sphere hit position and $\phi$
            pHit = ray((Float)tShapeHit);

            // Refine sphere intersection point
            pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
            if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
            phi = std::atan2(pHit.y, pHit.x);
            if (phi < 0) phi += 2 * Pi;
            if ((zMin > -radius && pHit.z < zMin) ||
                (zMax < radius && pHit.z > zMax) || phi > phiMax)
                return false;
        }
        return true;
    }

    Float Sphere::Area() const { return phiMax * radius * (zMax - zMin); }
    /*
    Interaction Sphere::Sample(const Point2f& u, Float* pdf) const {
        Point3f pObj = Point3f(0, 0, 0) + radius * UniformSampleSphere(u);
        Interaction it;
        it.n = Normalize((*ObjectToWorld)(Normal3f(pObj.x, pObj.y, pObj.z)));
        if (reverseOrientation) it.n *= -1;
        // Reproject _pObj_ to sphere surface and compute _pObjError_
        pObj *= radius / Distance(pObj, Point3f(0, 0, 0));
        Vector3f pObjError = gamma(5) * Abs((Vector3f)pObj);
        it.p = (*ObjectToWorld)(pObj, pObjError, &it.pError);
        *pdf = 1 / Area();
        return it;
    }

    Interaction Sphere::Sample(const Interaction& ref, const Point2f& u,
        Float* pdf) const {
        Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));

        // Sample uniformly on sphere if $\pt{}$ is inside it
        Point3f pOrigin =
            OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
        if (DistanceSquared(pOrigin, pCenter) <= radius * radius) {
            Interaction intr = Sample(u, pdf);
            Vector3f wi = intr.p - ref.p;
            if (wi.LengthSquared() == 0)
                *pdf = 0;
            else {
                // Convert from area measure returned by Sample() call above to
                // solid angle measure.
                wi = Normalize(wi);
                *pdf *= DistanceSquared(ref.p, intr.p) / AbsDot(intr.n, -wi);
            }
            if (std::isinf(*pdf)) *pdf = 0.f;
            return intr;
        }

        // Sample sphere uniformly inside subtended cone

        // Compute coordinate system for sphere sampling
        Float dc = Distance(ref.p, pCenter);
        Float invDc = 1 / dc;
        Vector3f wc = (pCenter - ref.p) * invDc;
        Vector3f wcX, wcY;
        CoordinateSystem(wc, &wcX, &wcY);

        // Compute $\theta$ and $\phi$ values for sample in cone
        Float sinThetaMax = radius * invDc;
        Float sinThetaMax2 = sinThetaMax * sinThetaMax;
        Float invSinThetaMax = 1 / sinThetaMax;
        Float cosThetaMax = std::sqrt(std::max((Float)0.f, 1 - sinThetaMax2));

        Float cosTheta = (cosThetaMax - 1) * u[0] + 1;
        Float sinTheta2 = 1 - cosTheta * cosTheta;

        if (sinThetaMax2 < 0.00068523f // sin^2(1.5 deg) //) {
            // Fall back to a Taylor series expansion for small angles, where
             //  the standard approach suffers from severe cancellation errors 
            sinTheta2 = sinThetaMax2 * u[0];
            cosTheta = std::sqrt(1 - sinTheta2);
        }

        // Compute angle $\alpha$ from center of sphere to sampled point on surface
        Float cosAlpha = sinTheta2 * invSinThetaMax +
            cosTheta * std::sqrt(std::max((Float)0.f, 1.f - sinTheta2 * invSinThetaMax * invSinThetaMax));
        Float sinAlpha = std::sqrt(std::max((Float)0.f, 1.f - cosAlpha * cosAlpha));
        Float phi = u[1] * 2 * Pi;

        // Compute surface normal and sampled point on sphere
        Vector3f nWorld =
            SphericalDirection(sinAlpha, cosAlpha, phi, -wcX, -wcY, -wc);
        Point3f pWorld = pCenter + radius * Point3f(nWorld.x, nWorld.y, nWorld.z);

        // Return _Interaction_ for sampled point on sphere
        Interaction it;
        it.p = pWorld;
        it.pError = gamma(5) * Abs((Vector3f)pWorld);
        it.n = Normal3f(nWorld);
        if (reverseOrientation) it.n *= -1;

        // Uniform cone PDF.
        *pdf = 1 / (2 * Pi * (1 - cosThetaMax));

        return it;
    }

    Float Sphere::Pdf(const Interaction& ref, const Vector3f& wi) const {
        Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));
        // Return uniform PDF if point is inside sphere
        Point3f pOrigin =
            OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
        if (DistanceSquared(pOrigin, pCenter) <= radius * radius)
            return Shape::Pdf(ref, wi);

        // Compute general sphere PDF
        Float sinThetaMax2 = radius * radius / DistanceSquared(ref.p, pCenter);
        Float cosThetaMax = std::sqrt(std::max((Float)0, 1 - sinThetaMax2));
        return UniformConePdf(cosThetaMax);
    }

    Float Sphere::SolidAngle(const Point3f& p, int nSamples) const {
        Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));
        if (DistanceSquared(p, pCenter) <= radius * radius)
            return 4 * Pi;
        Float sinTheta2 = radius * radius / DistanceSquared(p, pCenter);
        Float cosTheta = std::sqrt(std::max((Float)0, 1 - sinTheta2));
        return (2 * Pi * (1 - cosTheta));
    }

    std::shared_ptr<Shape> CreateSphereShape(const Transform* o2w,
        const Transform* w2o,
        bool reverseOrientation,
        const ParamSet& params) {
        Float radius = params.FindOneFloat("radius", 1.f);
        Float zmin = params.FindOneFloat("zmin", -radius);
        Float zmax = params.FindOneFloat("zmax", radius);
        Float phimax = params.FindOneFloat("phimax", 360.f);
        return std::make_shared<Sphere>(o2w, w2o, reverseOrientation, radius, zmin,
            zmax, phimax);
    }
    */
}  // namespace pbrt

