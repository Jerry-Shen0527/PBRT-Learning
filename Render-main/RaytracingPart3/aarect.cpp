#include "aarect.h"
bool xy_rect::Intersect(const ray& r, double t_min, double t_max, SurfaceInteraction& inter_record) const {
    auto t = (k - r.origin().z) / r.direction().z;
    if (t < t_min || t > t_max)
        return false;
    auto x = r.origin().x + t * r.direction().x;
    auto y = r.origin().y + t * r.direction().y;
    if (x < x0 || x > x1 || y < y0 || y > y1)
        return false;
    inter_record.u = (x - x0) / (x1 - x0);
    inter_record.v = (y - y0) / (y1 - y0);
    inter_record.t = t;
    auto outward_normal = Vector3f(0, 0, 1);
    inter_record.set_face_normal(r, outward_normal);
    inter_record.mat_ptr = mp;
    inter_record.p = r.at(t);
    return true;
}

bool xy_rect::IntersectP(const ray& r, double t_min, double t_max) const {
    auto t = (k - r.origin().z) / r.direction().z;
    if (t < t_min || t > t_max)
        return false;
    auto x = r.origin().x + t * r.direction().x;
    auto y = r.origin().y + t * r.direction().y;
    if (x < x0 || x > x1 || y < y0 || y > y1)
        return false;
    return true;
}

bool xz_rect::Intersect(const ray& r, double t_min, double t_max, SurfaceInteraction& inter_record) const {
    auto t = (k - r.origin()[1]) / r.direction()[1];
    //std::cout << t << " 1" << std::endl;
    //std::cout << r.origin()[1] << "  "<<k <<" "<< r.direction()[1]<<" "<<t <<std::endl;
    if (t < t_min || t > t_max)
        return false;
    auto x = r.origin().x + t * r.direction().x;
    auto z = r.origin().z + t * r.direction().z;

    
    //cout <<"111   "<< x << " " << t * r.direction().x << endl;
    //cout << r.origin().x << " " << r.origin().y << " " << r.origin().z << " " << r.direction().x << endl;
    //cout << r.direction().x << " " << r.direction().y << " " << r.direction().z << " " << r.direction().x << endl;
    if (x < x0 || x > x1 || z < z0 || z > z1)
        return false;
    //std::cout << x << " " << z << " " << x0 << " " << x1 << " " << z0 << " " << z1 << std::endl;
    //std::cout << r.origin().z << " " << t << " " << r.direction().z << std::endl;
    inter_record.u = (x - x0) / (x1 - x0);
    inter_record.v = (z - z0) / (z1 - z0);
    inter_record.t = t;
    //std::cout << t << " 2" << std::endl;
    //cout << x << " " << r.origin().x << " " << t * r.direction().x << " " << t << endl;
    auto outward_normal = Vector3f(0, 1, 0);
    inter_record.set_face_normal(r, outward_normal);
    inter_record.mat_ptr = mp;
    inter_record.p = r.at(t);

    //cout << rec.p.x << " " << rec.p.y << " " << rec.p.z << endl;
    return true;
}
bool xz_rect::IntersectP(const ray& r, double t_min, double t_max) const {
    auto t = (k - r.origin()[1]) / r.direction()[1];
    //std::cout << t << " 1" << std::endl;
    //std::cout << r.origin()[1] << "  "<<k <<" "<< r.direction()[1]<<" "<<t <<std::endl;
    if (t < t_min || t > t_max)
        return false;
    auto x = r.origin().x + t * r.direction().x;
    auto z = r.origin().z + t * r.direction().z;


    //cout <<"111   "<< x << " " << t * r.direction().x << endl;
    //cout << r.origin().x << " " << r.origin().y << " " << r.origin().z << " " << r.direction().x << endl;
    //cout << r.direction().x << " " << r.direction().y << " " << r.direction().z << " " << r.direction().x << endl;
    if (x < x0 || x > x1 || z < z0 || z > z1)
        return false;
    
    return true;
}

bool yz_rect::Intersect(const ray& r, double t_min, double t_max, SurfaceInteraction& inter_record) const {
    auto t = (k - r.origin().x) / r.direction().x;
    if (t < t_min || t > t_max)
        return false;
    auto y = r.origin().y + t * r.direction().y;
    auto z = r.origin().z + t * r.direction().z;
    if (y < y0 || y > y1 || z < z0 || z > z1)
        return false;
    inter_record.u = (y - y0) / (y1 - y0);
    inter_record.v = (z - z0) / (z1 - z0);
    inter_record.t = t;
    auto outward_normal = Vector3f(1, 0, 0);
    inter_record.set_face_normal(r, outward_normal);
    inter_record.mat_ptr = mp;
    inter_record.p = r.at(t);
    return true;
}
bool yz_rect::IntersectP(const ray& r, double t_min, double t_max) const {
    auto t = (k - r.origin().x) / r.direction().x;
    if (t < t_min || t > t_max)
        return false;
    auto y = r.origin().y + t * r.direction().y;
    auto z = r.origin().z + t * r.direction().z;
    if (y < y0 || y > y1 || z < z0 || z > z1)
        return false;

    return true;
}