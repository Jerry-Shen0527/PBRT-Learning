#ifndef MSPHERE_H
#define MSPHERE_H

#include "hittable.h"
#include "vec3.h"
#include "pdf.h"
#include "mtransform.h"

class msphere : public hittable {
public:
    msphere() {}
    msphere(Transform otw, Transform wto, double r, shared_ptr<material> m)
        : ObjectToWorld(otw), WorldToObject(wto), radius(r), mat_ptr(m) {};
    // const Transform*...

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override;

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override;

    double pdf_value(const point3& o, const vec3& v) const;

    vec3 random(const point3& o) const;
    //it seems no apply
public:
    point3 center;
    double radius;
    shared_ptr<material> mat_ptr;
    Transform ObjectToWorld, WorldToObject;

private:
    static void get_sphere_uv(const point3& p, double& u, double& v) {
        // p: a given point on the sphere of radius one, centered at the origin.
        // u: returned value [0,1] of angle around the Y axis from X=-1.
        // v: returned value [0,1] of angle from Y=-1 to Y=+1.
        //     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
        //     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
        //     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>

        auto theta = acos(-p.y());
        auto phi = atan2(-p.z(), p.x()) + pi;

        u = phi / (2 * pi);
        v = theta / pi;
    }
};

bool msphere::hit(const ray& wr, double t_min, double t_max, hit_record& rec) const {

    ray r = WorldToObject.ray_transform(wr);


    vec3 oc = r.origin() - center;

    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius * radius;

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
    
    r = ObjectToWorld.ray_transform(r);

    rec.t = root;
    rec.p = wr.at(rec.t);
    vec3 outward_normal = (rec.p - ObjectToWorld.point_transform(center)) / radius;
    rec.set_face_normal(wr, outward_normal);
    get_sphere_uv(outward_normal, rec.u, rec.v);
    rec.mat_ptr = mat_ptr;



    return true;
}

bool msphere::bounding_box(double time0, double time1, aabb& output_box) const {
    output_box = aabb(
        ObjectToWorld.point_transform(center) - vec3(radius, radius, radius),
        ObjectToWorld.point_transform(center) + vec3(radius, radius, radius));
    return true;
}

double msphere::pdf_value(const point3& o, const vec3& v) const {
    hit_record rec;
    if (!this->hit(ray(o, v), 0.001, infinity, rec))
        return 0;

    auto cos_theta_max = sqrt(1 - radius * radius / (center - o).length_squared());
    auto solid_angle = 2 * pi * (1 - cos_theta_max);

    return  1 / solid_angle;
}

vec3 msphere::random(const point3& o) const {
    vec3 direction = center - o;
    auto distance_squared = direction.length_squared();
    onb uvw;
    uvw.build_from_w(direction);
    return uvw.local(random_to_sphere(radius, distance_squared));
}
#endif
