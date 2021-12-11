#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "shape.h"

#include <map>

// Triangle Declarations
/*输入为局部坐标，而类的p n 等中存储世界坐标*/
struct TriangleMesh {
    // TriangleMesh Public Methods
    TriangleMesh(shared_ptr<Transform> ObjectToWorld, int nTriangles,
        const int* vertexIndices, int nVertices, const Point3f* P,
        const Vector3f* S, const Normal3f* N, const Point2f* uv,
        //const std::shared_ptr<Texture<Float>>& alphaMask,
        //const std::shared_ptr<Texture<Float>>& shadowAlphaMask,
        const int* faceIndices);

    // TriangleMesh Data
    const int nTriangles, nVertices;
    std::vector<int> vertexIndices;
    std::unique_ptr<Point3f[]> p;
    std::unique_ptr<Normal3f[]> n;
    std::unique_ptr<Vector3f[]> s;
    std::unique_ptr<Point2f[]> uv;
    //std::shared_ptr<Texture<Float>> alphaMask, shadowAlphaMask;
    std::vector<int> faceIndices;
};

/*存储mesh中第triNumber个三角形*/
class Triangle :public Shape {
public:
    // Triangle Public Methods
    Triangle(shared_ptr<Transform> ObjectToWorld, shared_ptr<Transform> WorldToObject,
        bool reverseOrientation, const std::shared_ptr<TriangleMesh>& mesh,
        int triNumber)
        : Shape(ObjectToWorld, WorldToObject, reverseOrientation), mesh(mesh) {
        v = &mesh->vertexIndices[3 * triNumber];
        //triMeshBytes += sizeof(*this);
        //faceIndex = mesh->faceIndices.size() ? mesh->faceIndices[triNumber] : 0;
    }

    Bounds3f ObjectBound() const;
    Bounds3f WorldBound() const;
    bool Intersect(const Ray& ray, Float* tHit, SurfaceInteraction* isect,
        bool testAlphaTexture = true) const;
    bool IntersectP(const Ray& ray, bool testAlphaTexture = true) const;
    Float Area() const;

    //using Shape::Sample;  // Bring in the other Sample() overload.
   // Interaction Sample(const Point2f& u, Float* pdf) const;

    // Returns the solid angle subtended by the triangle w.r.t. the given
    // reference point p.
    Float SolidAngle(const Point3f& p, int nSamples = 0) const;

private:
    // Triangle Private Methods
    void GetUVs(Point2f uv[3]) const {
        if (mesh->uv) {
            uv[0] = mesh->uv[v[0]];
            uv[1] = mesh->uv[v[1]];
            uv[2] = mesh->uv[v[2]];
        }
        else {
            uv[0] = Point2f(0, 0);
            uv[1] = Point2f(1, 0);
            uv[2] = Point2f(1, 1);
        }
    }

    // Triangle Private Data
    std::shared_ptr<TriangleMesh> mesh;
    const int* v;
    int faceIndex;
};

#endif // !TRIANGLE_H
