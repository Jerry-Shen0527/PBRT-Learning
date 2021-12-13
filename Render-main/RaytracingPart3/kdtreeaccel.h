#pragma once

#include "efloat.h"
#include "primitive_list.h"
#include "memory.h"
#include "check.h"

// KdTreeAccel Declarations
struct KdAccelNode;
struct BoundEdge;
class KdTreeAccel : public Primitive_list {
public:
	// KdTreeAccel Public Methods
	KdTreeAccel(std::vector<std::shared_ptr<Primitive>> p,
		int isectCost = 80, int traversalCost = 1,
		Float emptyBonus = 0.5, int maxPrims = 1, int maxDepth = -1);
	aabb WorldBound() const { return bounds; }
	~KdTreeAccel();
	bool Intersect(const ray& ray,  double t_min, double t_max, SurfaceInteraction& isect) const;
	bool IntersectP(const ray& ray, double t_min, double t_max) const;
	double pdf_value(const Point3f& o, const Vector3f& v) const;

private:
	// KdTreeAccel Private Methods
	void buildTree(int nodeNum, const aabb& bounds,
		const std::vector<aabb>& primBounds, int* primNums,
		int nprims, int depth,
		const std::unique_ptr<BoundEdge[]> edges[3], int* prims0,
		int* prims1, int badRefines = 0);

	// KdTreeAccel Private Data
	const int isectCost, traversalCost, maxPrims;
	const Float emptyBonus;
	std::vector<std::shared_ptr<Primitive>> primitives;
	std::vector<int> primitiveIndices;
	KdAccelNode* nodes;
	int nAllocedNodes, nextFreeNode;
	aabb bounds;
};

struct KdToDo {
	const KdAccelNode* node;
	Float tMin, tMax;
};

std::shared_ptr<KdTreeAccel> CreateKdTreeAccelerator(
	std::vector<std::shared_ptr<Primitive>> prims);
