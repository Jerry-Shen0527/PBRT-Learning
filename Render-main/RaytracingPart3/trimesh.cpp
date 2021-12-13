#include "trimesh.h"
aabb Trimesh::bounding_box() const {
	aabb result = tris[0]->bounding_box();
	for (int i = 1; i < tris.size(); i++) {
		result = Union(result, tris[i]->bounding_box());
	}
	return result;
}

bool Trimesh::Intersect(const ray& ray, double t_min, double t_max, SurfaceInteraction& inter_record) const {
	/*for (int i = 0; i < tris.size(); i++) {
		if (tris[i]->Intersect(ray, t_min, t_max, inter_record))
			return true;
	}
	return false;*/
	return accelerator->Intersect(ray, t_min, t_max, inter_record);
}

bool Trimesh::IntersectP(const ray& ray, double t_min, double t_max) const {
	/*for (int i = 0; i < tris.size(); i++) {
		if (tris[i]->IntersectP(ray, t_min, t_max))
			return true;
	}
	return false;*/
	return accelerator->IntersectP(ray, t_min, t_max);
}

double Trimesh::pdf_value(const Point3f& origin, const Vector3f& v) const {
	ray r_tmp(origin, v);
	SurfaceInteraction rec;
	/*for (int i = 0; i < tris.size(); i++) {
		if (tris[i]->Intersect(r_tmp, 0.001, infinity, rec))
			return tris[i]->pdf_value(origin, v);
	}*/
	return accelerator->pdf_value(origin, v);
}