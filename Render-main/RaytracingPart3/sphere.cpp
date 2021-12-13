#include "sphere.h"

bool sphere::Intersect(const ray& r, double t_min, double t_max, SurfaceInteraction& inter_record) const {
    //cout << center.x << endl;
    /*Vector3f oc = r.origin() - center;
    double phi;
    Point3f pHit,pRelative;
    // Transform _Ray_ to object space
    Vector3f oErr, dErr;
    //Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Compute quadratic sphere coefficients

    // Initialize _Edouble_ ray coordinate values
    EFloat ox(r.origin().x, oErr.x), oy(r.origin().y, oErr.y), oz(r.origin().z, oErr.z);
    EFloat dx(r.direction().x, dErr.x), dy(r.direction().y, dErr.y), dz(r.direction().z, dErr.z);
    EFloat a = dx * dx + dy * dy + dz * dz;
    EFloat b = 2 * (dx * ox + dy * oy + dz * oz)- 2 * (dx * center.x + dy * center.y + dz * center.z);
    EFloat c = (ox- center.x) * (ox - center.x) + (oy - center.y) * (oy - center.y) + (oz - center.z) * (oz - center.z) - EFloat(radius) * EFloat(radius);

    // Solve quadratic equation for _t_ values
    EFloat t0, t1;
    if (!Quadratic(a, b, c, &t0, &t1)) return false;

    // Check quadric shape _t0_ and _t1_ for nearest intersection
    if (t0.UpperBound() > t_max || t1.LowerBound() <= 0) return false;
    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() <= 0) {
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > t_max) return false;
    }

    // Compute sphere hit position and $\phi$


    pHit = r.at(tShapeHit.PreciseValue());
    pRelative = pHit - center;


    // Refine sphere intersection point
    pRelative *= radius / pRelative.Length();
    if (pRelative.x == 0 && pRelative.y == 0) pRelative.x = 1e-5f * radius;
    phi = std::atan2(pRelative.y, pRelative.x);
    if (phi < 0) phi += 2 * PI;

    double zMin = center.z - radius, zMax = center.z + radius;
    // Test sphere intersection against clipping parameters
    if ((zMin > -radius && pRelative.z < zMin) || (zMax < radius && pRelative.z > zMax) ||
        phi > phiMax) {
        if (tShapeHit == t1) return false;
        if (t1.UpperBound() > t_max) return false;
        tShapeHit = t1;
        // Compute sphere hit position and $\phi$
        pHit = r.at((Float)tShapeHit);

        // Refine sphere intersection point
        pHit *= radius / pRelative.Length();
        if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) phi += 2 * PI;
        if ((zMin > -radius && pHit.z < zMin) ||
            (zMax < radius && pHit.z > zMax) || phi > phiMax)
            return false;
    }

    // Find parametric representation of sphere hit
    double u = phi / phiMax;
    double theta = std::acos(check::Clamp(pRelative.z / radius, -1.0, 1.0));
    double v = (theta - thetaMin) / (thetaMax - thetaMin);

    // Compute sphere $\dpdu$ and $\dpdv$
    double zRadius = std::sqrt(pRelative.x * pRelative.x + pRelative.y * pRelative.y);
    double invZRadius = 1 / zRadius;
    double cosPhi = pRelative.x * invZRadius;
    double sinPhi = pRelative.y * invZRadius;
    Vector3f dpdu(-phiMax * pRelative.y, phiMax * pRelative.x, 0);
    Vector3f dpdv =
        (thetaMax - thetaMin) *
        Vector3f(pRelative.z * cosPhi, pRelative.z * sinPhi, -radius * std::sin(theta));

    // Compute sphere $\dndu$ and $\dndv$
    Vector3f d2Pduu = -phiMax * phiMax * Vector3f(pRelative.x, pRelative.y, 0);
    Vector3f d2Pduv =
        (thetaMax - thetaMin) * pRelative.z * phiMax * Vector3f(-sinPhi, cosPhi, 0.);
    Vector3f d2Pdvv = -(thetaMax - thetaMin) * (thetaMax - thetaMin) *
        Vector3f(pRelative.x, pRelative.y, pRelative.z);

    // Compute coefficients for fundamental forms
    double E = Dot(dpdu, dpdu);
    double F = Dot(dpdu, dpdv);
    double G = Dot(dpdv, dpdv);
    Vector3f N = Normalize(Cross(dpdu, dpdv));
    double e = Dot(N, d2Pduu);
    double f = Dot(N, d2Pduv);
    double g = Dot(N, d2Pdvv);


    Vector3f outward_normal(pRelative / radius);

    // Compute $\dndu$ and $\dndv$ from fundamental form coefficients
    double invEGF2 = 1 / (E * G - F * F);
    Vector3f dndu = Vector3f((f * F - e * G) * invEGF2 * dpdu +
        (e * F - f * E) * invEGF2 * dpdv);
    Vector3f dndv = Vector3f((g * F - f * G) * invEGF2 * dpdu +
        (f * F - g * E) * invEGF2 * dpdv);

    // Compute error bounds for sphere intersection
    Vector3f pError = gamma(5) * Abs((Vector3f)pRelative);

    inter_record.t = (double)tShapeHit;;
    inter_record.p = pRelative+center;
    //cout << center.x << " " << center.y << " " << center.z << " " << radius << endl;
    //cout << inter_record.p.y << endl;
    //Vector3f outward_normal(Vector3f(rec.p - center) / radius);
    inter_record.set_face_normal(r, outward_normal);
    //Point3f outward_normal_p = Point3f(rec.p - center) / radius;
    get_sphere_uv(outward_normal, inter_record.u, inter_record.v);
    inter_record.mat_ptr = mat_ptr;

    return true;*/


    Float phi;
    Point3f pHit;
    // Transform _Ray_ to object space
    Vector3f oErr, dErr;
    //ray ray = (*WorldToObject)(r, &oErr, &dErr);
    ray ray =ConvertRayTrans(WorldToObject,r, &oErr, &dErr, t_max, t_min);
    //ray ray(r.origin() - center, r.direction(), r.time());
    //ray.origin() -= center;
   // cout << ray.origin().x << endl;
    // Compute quadratic sphere coefficients

    // Initialize _EFloat_ ray coordinate values
    EFloat ox(ray.origin().x, oErr.x), oy(ray.origin().y, oErr.y), oz(ray.origin().z, oErr.z);
    EFloat dx(ray.direction().x, dErr.x), dy(ray.direction().y, dErr.y), dz(ray.direction().z, dErr.z);
    EFloat a = dx * dx + dy * dy + dz * dz;
    EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
    EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

    // Solve quadratic equation for _t_ values
    EFloat t0, t1;
    if (!Quadratic(a, b, c, &t0, &t1)) return false;

    // Check quadric shape _t0_ and _t1_ for nearest intersection
    if (t0.UpperBound() > t_max || t1.LowerBound() <= 0) return false;
    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() <= 0) {
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > t_max) return false;
    }

    pHit = ray.at(tShapeHit.PreciseValue());

    // Refine sphere intersection point
    pHit *= radius / (pHit-Point3f(0, 0, 0)).Length();
    if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
    phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0) phi += 2 * PI;

    // Test sphere intersection against clipping parameters
    if ((zMin > -radius && pHit.z < zMin) || (zMax < radius && pHit.z > zMax) ||
        phi > phiMax) {
        if (tShapeHit == t1) return false;
        if (t1.UpperBound() > t_max) return false;
        tShapeHit = t1;
        // Compute sphere hit position and $\phi$
        pHit = ray.at((Float)tShapeHit);

        // Refine sphere intersection point
        pHit *= radius / (pHit - Point3f(0, 0, 0)).Length();
        if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) phi += 2 * PI;
        if ((zMin > -radius && pHit.z < zMin) ||
            (zMax < radius && pHit.z > zMax) || phi > phiMax)
            return false;
    }

    // Find parametric representation of sphere hit
    Float u = phi / phiMax;
    Float theta = std::acos(check::Clamp(pHit.z / radius, -1, 1));
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
    double invEGF2 = 1 / (E * G - F * F);
    Vector3f dndu = Vector3f((f * F - e * G) * invEGF2 * dpdu +
        (e * F - f * E) * invEGF2 * dpdv);
    Vector3f dndv = Vector3f((g * F - f * G) * invEGF2 * dpdu +
        (f * F - g * E) * invEGF2 * dpdv);

    // Compute error bounds for sphere intersection
    Vector3f pError = gamma(5) * Abs((Vector3f)pHit);

    inter_record.t = (double)tShapeHit;;
    inter_record.p = ConvertPTrans(ObjectToWorld,pHit, (pError));
    //inter_record.p = pHit + center;
    //cout << pHit.x << " " << inter_record.p.x << endl;
    Vector3f outward_normal(Vector3f(pHit - Point3f(0, 0, 0)) / radius);
    Vector3f trans_N = ConvertVTrans(ObjectToWorld, outward_normal);
    inter_record.set_face_normal(r, trans_N);
    //Point3f outward_normal_p = Point3f(rec.p - center) / radius;
    get_sphere_uv((Point3f)trans_N, inter_record.u, inter_record.v);
    inter_record.mat_ptr = mat_ptr;
    //cout << "1" << endl;
    return true;
}


bool sphere::IntersectP(const ray& r, double t_min, double t_max) const {
    //center = ConvertPTrans(ObjectToWorld, Point3f(0, 0, 0));
    /*Vector3f oc = r.origin() - center;
    
    double phi;
    Point3f pHit, pRelative;
    // Transform _Ray_ to object space
    Vector3f oErr, dErr;
    //Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Compute quadratic sphere coefficients

    // Initialize _Edouble_ ray coordinate values
    EFloat ox(r.origin().x, oErr.x), oy(r.origin().y, oErr.y), oz(r.origin().z, oErr.z);
    EFloat dx(r.direction().x, dErr.x), dy(r.direction().y, dErr.y), dz(r.direction().z, dErr.z);
    EFloat a = dx * dx + dy * dy + dz * dz;
    EFloat b = 2 * (dx * ox + dy * oy + dz * oz) - 2 * (dx * center.x + dy * center.y + dz * center.z);
    EFloat c = (ox - center.x) * (ox - center.x) + (oy - center.y) * (oy - center.y) + (oz - center.z) * (oz - center.z) - EFloat(radius) * EFloat(radius);

    // Solve quadratic equation for _t_ values
    EFloat t0, t1;
    if (!Quadratic(a, b, c, &t0, &t1)) return false;

    // Check quadric shape _t0_ and _t1_ for nearest intersection
    if (t0.UpperBound() > t_max || t1.LowerBound() <= 0) return false;
    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() <= 0) {
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > t_max) return false;
    }

    // Compute sphere hit position and $\phi$


    pHit = r.at(tShapeHit.PreciseValue());
    pRelative = pHit - center;


    // Refine sphere intersection point
    pRelative *= radius / pRelative.Length();
    if (pRelative.x == 0 && pRelative.y == 0) pRelative.x = 1e-5f * radius;
    phi = std::atan2(pRelative.y, pRelative.x);
    if (phi < 0) phi += 2 * PI;

    double zMin = center.z - radius, zMax = center.z + radius;
    // Test sphere intersection against clipping parameters
    if ((zMin > -radius && pRelative.z < zMin) || (zMax < radius && pRelative.z > zMax) ||
        phi > phiMax) {
        if (tShapeHit == t1) return false;
        if (t1.UpperBound() > t_max) return false;
        tShapeHit = t1;
        // Compute sphere hit position and $\phi$
        pHit = r.at((Float)tShapeHit);

        // Refine sphere intersection point
        pHit *= radius / pRelative.Length();
        if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) phi += 2 * PI;
        if ((zMin > -radius && pHit.z < zMin) ||
            (zMax < radius && pHit.z > zMax) || phi > phiMax)
            return false;
    }

    return true;*/




    
    Float phi;
    Point3f pHit;
    // Transform _Ray_ to object space
    Vector3f oErr, dErr;
    //ray ray = (*WorldToObject)(r, &oErr, &dErr);
    ray ray = ConvertRayTrans(WorldToObject, r, &oErr, &dErr, t_max, t_min);

    // Compute quadratic sphere coefficients

    // Initialize _EFloat_ ray coordinate values
    EFloat ox(ray.origin().x, oErr.x), oy(ray.origin().y, oErr.y), oz(ray.origin().z, oErr.z);
    EFloat dx(ray.direction().x, dErr.x), dy(ray.direction().y, dErr.y), dz(ray.direction().z, dErr.z);
    EFloat a = dx * dx + dy * dy + dz * dz;
    EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
    EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

    // Solve quadratic equation for _t_ values
    EFloat t0, t1;
    if (!Quadratic(a, b, c, &t0, &t1)) return false;

    // Check quadric shape _t0_ and _t1_ for nearest intersection
    if (t0.UpperBound() > t_max || t1.LowerBound() <= 0) return false;
    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() <= 0) {
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > t_max) return false;
    }

    pHit = ray.at(tShapeHit.PreciseValue());

    // Refine sphere intersection point
    pHit *= radius / (pHit - Point3f(0, 0, 0)).Length();
    if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
    phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0) phi += 2 * PI;

    // Test sphere intersection against clipping parameters
    if ((zMin > -radius && pHit.z < zMin) || (zMax < radius && pHit.z > zMax) ||
        phi > phiMax) {
        if (tShapeHit == t1) return false;
        if (t1.UpperBound() > t_max) return false;
        tShapeHit = t1;
        // Compute sphere hit position and $\phi$
        pHit = ray.at((Float)tShapeHit);

        // Refine sphere intersection point
        pHit *= radius / (pHit - Point3f(0, 0, 0)).Length();
        if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) phi += 2 * PI;
        if ((zMin > -radius && pHit.z < zMin) ||
            (zMax < radius && pHit.z > zMax) || phi > phiMax)
            return false;
    }

    return true;
}

aabb sphere::bounding_box() const {
    //cout << (center - Vector3f(radius, radius, radius)).x << endl;
    return aabb(
        center - Vector3f(radius, radius, radius),
        center + Vector3f(radius, radius, radius));
    
}

void sphere::get_sphere_uv(const Point3f& p, double& u, double& v) {
    // p: a given point on the sphere of radius one, centered at the origin.
    // u: returned value [0,1] of angle around the Y axis from X=-1.
    // v: returned value [0,1] of angle from Y=-1 to Y=+1.
    //     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
    //     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
    //     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>

    auto theta = acos(-p.y);
    auto phi = atan2(-p.z, p.x) + PI;

    u = phi / (2 * PI);
    v = theta / PI;
}

double sphere::pdf_value(const Point3f& o, const Vector3f& v) const {
    SurfaceInteraction inter_record;
    if (!this->Intersect(ray(o, v), 0.001, infinity, inter_record))
        return 0;

    auto cos_theta_max = sqrt(1 - radius * radius / (Vector3f(center - o)).LengthSquared());
    auto solid_angle = 2 * PI * (1 - cos_theta_max);

    return  1 / solid_angle;
}

Vector3f sphere::random(const Point3f& o) const {
    Vector3f direction(center - o);
    auto distance_squared = direction.LengthSquared();
    onb uvw;
    uvw.build_from_w(direction);
    return uvw.local(random_to_sphere(radius, distance_squared));
}