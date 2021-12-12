#ifndef MTRIANGLE_H
#define MTRIANGLE_H

#include "hittable.h"
#include "vec3.h"
#include "pdf.h"
#include "mtransform.h"
#include <vector>
#include "hittable_list.h"
using std::vector;

class triangle :public hittable {
public:
    triangle() {}
    triangle(point3 u1, point3 u2, point3 u3, shared_ptr<material> m)
        : v0(u1), v1(u2), v2(u3), mat_ptr(m) {

        normal = cross((v1 - v0), (v2 - v0));
        tri_area = 0.5 * normal.length();
        normal = unit_vector(normal);
    };
    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override;

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override;

    double pdf_value(const point3& o, const vec3& v) const;

    vec3 random(const point3& o) const;
    //no apply 

public:
    point3 v0, v1, v2;
    shared_ptr<material> mat_ptr;
    vec3 normal;
    double tri_area;

};

bool triangle::hit(const ray& r, double t_min, double t_max, hit_record& rec) const
{
    //r.hit_tri.hitted_triangle = true;
    
    point3 p0t = v0 - r.orig;
    point3 p1t = v1 - r.orig;
    point3 p2t = v2 - r.orig;

    int kz = maxdimension(vector_abs(r.dir));
    int kx = kz + 1;
    if (kx == 3) kx = 0;
    int ky = kx + 1;
    if (ky == 3)ky = 0;
    vec3 d = permutes(r.dir, kx, ky, kz);
    double Sx = -d.e[0] / d.e[2];
    double Sy = -d.e[1] / d.e[2];
    double Sz = 1.f / d.e[2];
    
    p0t = permutes(p0t, kx, ky, kz);
    p1t = permutes(p1t, kx, ky, kz);
    p2t = permutes(p2t, kx, ky, kz);

    
    p0t.e[0] += Sx * p0t.e[2];
    p0t.e[1] += Sy * p0t.e[2];
    p1t.e[0] += Sx * p1t.e[2];
    p1t.e[1] += Sy * p1t.e[2];
    p2t.e[0] += Sx * p2t.e[2];
    p2t.e[1] += Sy * p2t.e[2];
    // 上面的计算可以放到ray类中

    double e0 = p1t.e[0] * p2t.e[1] - p1t.e[1] * p2t.e[0];
    double e1 = p2t.e[0] * p0t.e[1] - p2t.e[1] * p0t.e[0];
    double e2 = p0t.e[0] * p1t.e[1] - p0t.e[1] * p1t.e[0];

    if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
    {
        return false;
    }

    double det = e0 + e1 + e2;
    if (det == 0)
    {
        return false;
    }

    p0t.e[2] *= Sz;
    p1t.e[2] *= Sz;
    p2t.e[2] *= Sz;

    double tScaled = e0 * p0t.e[2] + e1 * p1t.e[2] + e2 * p2t.e[2];
    if (det < 0 && (tScaled > t_min * det || tScaled < t_max * det))
    {
        return false;
    }
    else if (det > 0 && (tScaled<t_min * det || tScaled>t_max * det))
    {
        return false;
    }
    double invDet = 1 / det;
    double b0 = e0 * invDet;
    double b1 = e1 * invDet;
    double b2 = e2 * invDet;
    double t = tScaled * invDet;

    rec.t = t;
    rec.p = r.at(t);
    //rec.p = b0 * v0 + b1 * v1 + b2 * v2;
    rec.set_face_normal(r, this->normal);
    //uv 以v0为(0,0) v0v1-u v0v2-v
    rec.u = b1;
    rec.v = b2;
    rec.mat_ptr = mat_ptr;

    return true;


}

double triangle::pdf_value(const point3& origin, const vec3& v) const
{
    hit_record rec;
    if (!this->hit(ray(origin, v), 0.001, infinity, rec))
        return 0;

    auto distance_squared = rec.t * rec.t * v.length_squared();
    auto cosine = fabs(dot(v, rec.normal) / v.length());

    return distance_squared / (cosine * tri_area);
}

bool triangle::bounding_box(double time0, double time1, aabb& output_box) const
{
    point3 pmin = min(min(v0, v1), v1);
    point3 pmax = max(max(v0, v1), v2);
    output_box = aabb(pmin - 0.001 * vec3(1, 1, 1), pmax + 0.001 * vec3(1, 1, 1));
    return true;

    

}

vec3 triangle::random(const point3& origin) const
{
    auto random_point = v0+random_double(0,1)*(v1-v0)+random_double(0,1)*(v2-v0);
    return random_point - origin;
}


class trianglemesh :public hittable {
public:
    trianglemesh(){}
    trianglemesh(Transform otw, Transform wto, vector<point3> varray, vector<vector<int>> farray
        , shared_ptr<material> m)
    {
        ObjectToWorld = otw;
        WorldToObject = wto;
        Varray.assign(varray.begin(), varray.end());
        Farray.assign(farray.begin(), farray.end());
        mat_ptr = m;
        for (int i = 0; i < Farray.size(); i++)
        {
            //shared_ptr<hittable> triangle3 = make_shared<triangle>(v01, v11, v21, green);
            mesh.add(make_shared<triangle>(Varray[Farray[i][0]], Varray[Farray[i][1]], Varray[Farray[i][2]], mat_ptr));
        }
    }

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override;
    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override;
    hittable_list mesh_obj() { return mesh; }

public:
    //vector<int> verticesIndex;
    vector<point3> Varray;
    vector<vector<int>> Farray;
    Transform ObjectToWorld, WorldToObject;
    hittable_list mesh;

    shared_ptr<material> mat_ptr;
};

bool trianglemesh::hit(const ray& r, double t_min, double t_max, hit_record& rec) const
{
    ray Oray = WorldToObject.ray_transform(r);
    bool hitornot = mesh.hit(Oray, t_min, t_max, rec);

    if (hitornot)
    {
        rec.p = r.at(rec.t);
        rec.normal = ObjectToWorld.vector_transform(rec.normal);
        return true;
    }
    return false;
    
}

bool trianglemesh::bounding_box(double time0, double time1, aabb& output_box) const
{
    point3 small=Varray[0];
    point3 big=Varray[0];
    for (int i = 0; i < Varray.size(); i++)
    {
        small = min(small, Varray[i]);
        big = max(big, Varray[i]);
    }
    output_box = aabb(ObjectToWorld, WorldToObject,small, big);
    return true;
}



#endif