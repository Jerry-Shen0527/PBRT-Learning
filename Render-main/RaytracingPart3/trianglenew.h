#pragma once
#include "rtweekend.h"

#include "primitive.h"
#include "transform.h"
#define USE_OPENMESH
#ifdef USE_OPENMESH
#include<OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include<OpenMesh/Core/Mesh/Handles.hh>
#include<OpenMesh/Core/Mesh/IteratorsT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
typedef OpenMesh::PolyMesh_ArrayKernelT<> Mesh;
#endif // USE_OPENMESH

using namespace std;
class Trianglenew:public Primitive
{
public:
    Trianglenew(const Transform& OToWorld, const Transform& WToObject,
        bool reverse, const Point3f* pos,shared_ptr<material> m)
        : ObjectToWorld(OToWorld),
        WorldToObject(WToObject),
        reverseOrientation(reverse),
        mat_ptr(m){
        p[0] = pos[0];
        p[1] = pos[1];
        p[2] = pos[2];
        //v = &mesh->vertexIndices[3 * triNumber];
        ////triMeshBytes += sizeof(*this);
        //faceIndex = mesh->faceIndices.size() ? mesh->faceIndices[triNumber] : 0;
    }
    //aabb ObjectBound() const;
    virtual aabb bounding_box() const override;

    virtual bool Intersect(const ray& ray, double t_min, double t_max, SurfaceInteraction& inter_record) const override;
    virtual bool IntersectP(const ray& ray, double t_min, double t_max) const override;

    virtual double pdf_value(const Point3f& origin, const Vector3f& v) const override {
        SurfaceInteraction rec;
        if (!this->Intersect(ray(origin, v), 0.001, infinity, rec))
            return 0;

        auto area = Area();
        //cout << rec.t << "1" << v.LengthSquared() << endl;
        auto distance_squared = rec.t * rec.t * v.LengthSquared();
        //std::cout << "2 " << rec.t << std::endl;
        auto cosine = fabs(Dot(v, rec.normal) / v.Length());
        //std::cout << distance_squared<<" "<<(cosine)<<" "<<  area << std::endl;
        return distance_squared / (cosine * area);
    }

    virtual Vector3f random(const Point3f& origin) const override {
        double r1 = random_double(0.001, 0.999);
        double r2 = random_double(0.001, 0.999);
        double r3 = random_double(0.001, 0.999);
        double rs = r1 + r2 + r3;
        r1 /= rs;
        r2 /= rs;
        r3 /= rs;
        auto random_point = r1 * p[0] + r2 * p[1] + r3 * p[2];
        return Vector3f(random_point - origin);
    }

    Float Area() const {
        /*const Point3f& p0 = mesh->p[v[0]];
        const Point3f& p1 = mesh->p[v[1]];
        const Point3f& p2 = mesh->p[v[2]];*/
        return 0.5 * Cross(p[1] - p[0], p[2] - p[0]).Length();
    }
private:   
    void GetUVs(Point2f uv[3]) const {
        /*if (mesh->uv) {
            uv[0] = mesh->uv[v[0]];
            uv[1] = mesh->uv[v[1]];
            uv[2] = mesh->uv[v[2]];
        }
        else {*/
            uv[0] = Point2f(0, 0);
            uv[1] = Point2f(1, 0);
            uv[2] = Point2f(1, 1);
        //}
    }
    // Triangle Private Data
    const Transform ObjectToWorld, WorldToObject;
    const bool reverseOrientation;
    Point3f p[3];
    //const int* v;
    //int faceIndex;
    shared_ptr<material> mat_ptr;
};

struct TriangleMeshNew {
    // TriangleMesh Public Methods
    TriangleMeshNew(const Transform& ObjectToWorld, bool reverseOrientation,const int nTriangles,
        const int nVertices, const vector<Point3f>& v, vector<int>& f, vector< shared_ptr<material>>& m);

    // TriangleMesh Data
    const int nTriangles, nVertices;
    //std::vector<int> vertexIndices;
    //std::unique_ptr<Point3f[]> p;
    //std::unique_ptr<Vector3f[]> n;
    //std::unique_ptr<Vector3f[]> s;
    //std::unique_ptr<Point2f[]> uv;
    ////std::shared_ptr<Texture<Float>> alphaMask, shadowAlphaMask;
    //std::vector<int> faceIndices;
    const Transform ObjectToWorld;
    const bool reverseOrientation;
    vector<int> faceVID;//for each face, the vertex id of the face is faceVID[3*i], faceVID[3*i+1], faceVID[3*i+2]
    vector<Point3f> vPos;
    vector< shared_ptr<material>> mat_ptr;
};

std::vector<std::shared_ptr<Primitive>> CreateTriangleMesh(
    const Transform& o2w, const Transform& w2o, bool reverseOrientation,
    const int nTriangles,const int nVertices,const vector<Point3f>& vPos, vector<int>& faceVID, vector< shared_ptr<material>>& m);

std::vector<std::shared_ptr<Primitive>> CreateTriangleMesh(
    const Transform& o2w, const Transform& w2o, bool reverseOrientation,
    const int nTriangles, const int nVertices, const vector<Point3f>& vPos, vector<int>& faceVID, shared_ptr<material> m);

void ReadObj(const char* filename,vector<Point3f>& vPos, vector<int>& faceVID,double scale);

#ifdef USE_OPENMESH
void ReadObjOpenMesh(const char* filename, vector<Point3f>& vPos, vector<int>& faceVID, double scale);
#endif // USE_OPENMESH
