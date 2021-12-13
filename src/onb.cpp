#include "onb.h"


void onb::build_from_w(const Normal3f& n) {
    axis[2] = unit_vector(Vector3f(n));
    Vector3f a = (fabs(w().x) > 0.9) ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0);
    axis[1] = unit_vector(Cross(w(), a));
    axis[0] = Cross(w(), v());
}

void onb::build_from_w(const Vector3f& n) {
    axis[2] = unit_vector(n);
    Vector3f a = (fabs(w().x) > 0.9) ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0);
    axis[1] = unit_vector(Cross(w(), a));
    axis[0] = Cross(w(), v());
}