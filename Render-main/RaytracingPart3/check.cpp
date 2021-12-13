#include "check.h"

void check::CHECK(bool event) {
    if (!event) std::cout << "Check failed!" << std::endl;
}

void check::DCHECK(bool event) {
    if (!event) std::cout << "Check failed!" << std::endl;
}

void check::CHECK_LT(double v1, double v2) {
    if (!(v1 < v2)) std::cout << "Check failed!: CHECK_LT" << std::endl;
}

void check::CHECK_GT(double v1, double v2) {
    if (!(v1 > v2)) std::cout << "Check failed!: CHECK_GT" << std::endl;
}

void check::CHECK_GE(double v1, double v2) {
    if (!(v1 >= v2)) std::cout << "Check failed!: CHECK_GE" << std::endl;
}

void check::CHECK_LE(double v1, double v2) {
    if (!(v1 <= v2)) std::cout << "Check failed!: CHECK_LE" << std::endl;
}

void check::CHECK_EQ(double v1, double v2) {
    if (!(v1 == v2)) std::cout << "Check failed!: CHECK_EQ" << std::endl;
}

void check::CHECK_NE(double v1, double v2) {
    if (!(v1 != v2)) std::cout << "Check failed!: CHECK_NE" << std::endl;
}


template <typename T, typename U, typename V>
T check::Clamp(T val, U low, V high) {
    if (val < low)
        return low;
    else if (val > high)
        return high;
    else
        return val;
}

template int check::Clamp<int>(int, int, int);
template double check::Clamp<double>(double, double, double);
template double check::Clamp<double>(double, int, int);
template float check::Clamp<float>(float, float, float);
template float check::Clamp<float>(float, int, int);

double check::Lerp(double t, double v1, double v2) { return (1 - t) * v1 + t * v2; }