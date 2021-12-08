#include "hittable.h"
#include "shape.h"
#include "spectrum.h"

bool translate::hit(const ray& r, Float t_min, Float t_max, hit_record& rec) const {
    ray moved_r(r.origin() - offset, r.direction(), r.Time());
    if (!ptr->hit(moved_r, t_min, t_max, rec))
        return false;

    rec.p += offset;
    rec.set_face_normal(moved_r, rec.normal);

    return true;
}

bool translate::bounding_box(Float time0, Float time1, aabb& output_box) const {
    if (!ptr->bounding_box(time0, time1, output_box))
        return false;

    output_box = aabb(
        output_box.min() + offset,
        output_box.max() + offset);

    return true;
}

 bool rotate_y::bounding_box(Float time0, Float time1, aabb& output_box) const {
    output_box = bbox;
    return hasbox;
}

rotate_y::rotate_y(shared_ptr<hittable> p, Float angle) : ptr(p) {
    auto radians = Radians(angle);
    sin_theta = sin(radians);
    cos_theta = cos(radians);
    hasbox = ptr->bounding_box(0, 1, bbox);

    point3 min(Infinity, Infinity, Infinity);
    point3 max(-Infinity, -Infinity, -Infinity);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                auto x = i * bbox.max().x + (1 - i) * bbox.min().x;
                auto y = j * bbox.max().y + (1 - j) * bbox.min().y;
                auto z = k * bbox.max().z + (1 - k) * bbox.min().z;

                auto newx = cos_theta * x + sin_theta * z;
                auto newz = -sin_theta * x + cos_theta * z;

                vec3 tester(newx, y, newz);

                for (int c = 0; c < 3; c++) {
                    min[c] = fmin(min[c], tester[c]);
                    max[c] = fmax(max[c], tester[c]);
                }
            }
        }
    }

    bbox = aabb(min, max);
}

bool rotate_y::hit(const ray& r, Float t_min, Float t_max, hit_record& rec) const {
    auto origin = r.origin();
    auto direction = r.direction();

    origin[0] = cos_theta * r.origin()[0] - sin_theta * r.origin()[2];
    origin[2] = sin_theta * r.origin()[0] + cos_theta * r.origin()[2];

    direction[0] = cos_theta * r.direction()[0] - sin_theta * r.direction()[2];
    direction[2] = sin_theta * r.direction()[0] + cos_theta * r.direction()[2];

    ray rotated_r(origin, direction, r.Time());

    if (!ptr->hit(rotated_r, t_min, t_max, rec))
        return false;

    auto p = rec.p;
    auto normal = rec.normal;

    p[0] = cos_theta * rec.p[0] + sin_theta * rec.p[2];
    p[2] = -sin_theta * rec.p[0] + cos_theta * rec.p[2];

    normal[0] = cos_theta * rec.normal[0] + sin_theta * rec.normal[2];
    normal[2] = -sin_theta * rec.normal[0] + cos_theta * rec.normal[2];

    rec.p = p;
    rec.set_face_normal(rotated_r, normal);

    return true;
}



 SurfaceInteraction::SurfaceInteraction(const point3& p, const vec3& pError, Point2f uv, const vec3& wo, const vec3& dpdu, const vec3& dpdv, const Normal& dndu, const Normal& dndv, Float time, const Shape* shape, int faceIndex)
    : hit_record(p, Normal(Normalize(Cross(dpdu, dpdv))), pError, wo, time),
    uv(uv),
    dpdu(dpdu),
    dpdv(dpdv),
    dndu(dndu),
    dndv(dndv),
    shape(shape)
{
    // Initialize shading geometry from true geometry
    shading.n = n;
    shading.dpdu = dpdu;
    shading.dpdv = dpdv;
    shading.dndu = dndu;
    shading.dndv = dndv;

    // Adjust normal based on orientation and handedness
    if (shape &&
        (shape->reverseOrientation ^ shape->transformSwapsHandedness)) {
        n *= -1;
        shading.n *= -1;
    }
}

 void SurfaceInteraction::SetShadingGeometry(const vec3& dpdus,
     const vec3& dpdvs, const Normal& dndus,
     const Normal& dndvs, bool orientationIsAuthoritative) {
     shading.n = Normalize((Normal)Cross(dpdus, dpdvs));
     if (shape && (shape->reverseOrientation ^
         shape->transformSwapsHandedness))
         shading.n = -shading.n;
     if (orientationIsAuthoritative)
         n = Faceforward(n, vec3(shading.n));
     else
         shading.n = Faceforward(shading.n, vec3(n));

     shading.dpdu = dpdus;
     shading.dpdv = dpdvs;
     shading.dndu = dndus;
     shading.dndv = dndvs;

 }

 RGBSpectrum SurfaceInteraction::Le(const Vector3f& w) const {
     //const AreaLight* area = primitive->GetAreaLight();
     //return area ? area->L(*this, w) : RGBSpectrum(0.f);
     return RGBSpectrum(0.f);
 }
