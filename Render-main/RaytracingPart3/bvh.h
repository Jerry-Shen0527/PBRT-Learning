#pragma once
#include "primitive_list.h"
#include "memory.h"
struct BVHBuildNode;

// BVHAccel Forward Declarations
struct BVHPrimitiveInfo;
struct MortonPrimitive;
struct LinearBVHNode;

// BVHAccel Declarations
class BVHAccel : public Primitive_list {
public:
    // BVHAccel Public Types
    enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts };

    // BVHAccel Public Methods
    BVHAccel(std::vector<std::shared_ptr<Primitive>> p,
        int maxPrimsInNode = 1,
        SplitMethod splitMethod = SplitMethod::SAH);
    virtual aabb bounding_box() const override; 
    ~BVHAccel();
    bool Intersect(const ray& ray, double t_min, double t_max, SurfaceInteraction& isect) const;
    bool IntersectP(const ray& ray, double t_min, double t_max) const;
    double pdf_value(const Point3f& o, const Vector3f& v) const;

private:
    // BVHAccel Private Methods
    BVHBuildNode* recursiveBuild(
        MemoryArena& arena, std::vector<BVHPrimitiveInfo>& primitiveInfo,
        int start, int end, int* totalNodes,
        std::vector<std::shared_ptr<Primitive>>& orderedPrims);
    BVHBuildNode* HLBVHBuild(
        MemoryArena& arena, const std::vector<BVHPrimitiveInfo>& primitiveInfo,
        int* totalNodes,
        std::vector<std::shared_ptr<Primitive>>& orderedPrims) const;
    BVHBuildNode* emitLBVH(
        BVHBuildNode*& buildNodes,
        const std::vector<BVHPrimitiveInfo>& primitiveInfo,
        MortonPrimitive* mortonPrims, int nPrimitives, int* totalNodes,
        std::vector<std::shared_ptr<Primitive>>& orderedPrims,
        std::atomic<int>* orderedPrimsOffset, int bitIndex) const;
    BVHBuildNode* buildUpperSAH(MemoryArena& arena,
        std::vector<BVHBuildNode*>& treeletRoots,
        int start, int end, int* totalNodes) const;
    int flattenBVHTree(BVHBuildNode* node, int* offset);

    // BVHAccel Private Data
    const int maxPrimsInNode;
    const SplitMethod splitMethod;
    std::vector<std::shared_ptr<Primitive>> primitives;
    LinearBVHNode* nodes = nullptr;
};

std::shared_ptr<BVHAccel> CreateBVHAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims);
