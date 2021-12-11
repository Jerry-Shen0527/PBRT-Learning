#include "bvh.h"

BVHAccel::BVHAccel(std::vector<std::shared_ptr<Primitive>> p,
    int maxPrimsInNode, SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)),
    splitMethod(splitMethod),
    primitives(std::move(p)) {
    //ProfilePhase _(Prof::AccelConstruction);
    if (primitives.empty()) return;
    // Build BVH from _primitives_

    //// Initialize _primitiveInfo_ array for primitives
    //std::vector<BVHPrimitiveInfo> primitiveInfo(primitives.size());
    //for (size_t i = 0; i < primitives.size(); ++i)
    //    primitiveInfo[i] = { i, primitives[i]->WorldBound() };

    //// Build BVH tree for primitives using _primitiveInfo_
    ////MemoryArena arena(1024 * 1024);
    //int totalNodes = 0;
    //std::vector<std::shared_ptr<Primitive>> orderedPrims;
    //orderedPrims.reserve(primitives.size());
    //BVHBuildNode* root;
    //if (splitMethod == SplitMethod::HLBVH)
    //    root = HLBVHBuild(arena, primitiveInfo, &totalNodes, orderedPrims);
    //else
    //    root = recursiveBuild(arena, primitiveInfo, 0, primitives.size(),
    //        &totalNodes, orderedPrims);
    //primitives.swap(orderedPrims);
    //primitiveInfo.resize(0);
    //LOG(INFO) << StringPrintf("BVH created with %d nodes for %d "
    //    "primitives (%.2f MB), arena allocated %.2f MB",
    //    totalNodes, (int)primitives.size(),
    //    float(totalNodes * sizeof(LinearBVHNode)) /
    //    (1024.f * 1024.f),
    //    float(arena.TotalAllocated()) /
    //    (1024.f * 1024.f));

    //// Compute representation of depth-first traversal of BVH tree
    //treeBytes += totalNodes * sizeof(LinearBVHNode) + sizeof(*this) +
    //    primitives.size() * sizeof(primitives[0]);
    //nodes = AllocAligned<LinearBVHNode>(totalNodes);
    //int offset = 0;
    //flattenBVHTree(root, &offset);
    //CHECK_EQ(totalNodes, offset);
}
