#pragma once

#include "trianglenew.h"
#include "bvh.h"
using namespace std;

class Trimesh:public Primitive
{
public:
    Trimesh(const Transform& OToWorld, const Transform& WToObject,
        bool reverse, const vector<Point3f>& pos, const vector<int>& face, shared_ptr<material> m)
        : ObjectToWorld(OToWorld),
        WorldToObject(WToObject),
        reverseOrientation(reverse),
        vPos(pos),
        faceVID(face),
        mat_ptr(m) {
        for (int i = 0; i < face.size() / 3; i++) {
            Point3f p[3] = { ConvertPTrans(OToWorld,pos[face[3 * i]]),ConvertPTrans(OToWorld,pos[face[3 * i + 1]]) 
                ,ConvertPTrans(OToWorld,pos[face[3 * i + 2]]) };
            tris.push_back(make_shared<Trianglenew>(OToWorld, WToObject, reverse, p, m));
        }
        accelerator =CreateBVHAccelerator(tris);
    }
    //aabb ObjectBound() const;
    virtual aabb bounding_box() const override;

    virtual bool Intersect(const ray& ray, double t_min, double t_max, SurfaceInteraction& inter_record) const override;
    virtual bool IntersectP(const ray& ray, double t_min, double t_max) const override;

    virtual double pdf_value(const Point3f& origin, const Vector3f& v) const override;

private:
    const Transform ObjectToWorld, WorldToObject;
    const bool reverseOrientation;
    vector<Point3f> vPos;
    vector<int> faceVID;
    shared_ptr<material> mat_ptr;
    vector<shared_ptr<Primitive>> tris;
    shared_ptr<Primitive> accelerator;
};

