#ifndef PBRT_H
#define PBRT_H

#include <type_traits>
#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <assert.h>
#include <string.h>

#include <stdint.h>

#include <float.h>
#include <intrin.h>



namespace pbrt {

    template <typename T>
    class Vector2;
    template <typename T>
    class Vector3;
    template <typename T>
    class Point3;
    template <typename T>
    class Point2;
    template <typename T>
    class Normal3;

    class RGBSpectrum;
    template <int nSpectrumSamples>
    class CoefficientSpectrum;
    class SampledSpectrum;
#ifdef PBRT_SAMPLED_SPECTRUM
    typedef SampledSpectrum Spectrum;
#else
    typedef RGBSpectrum Spectrum;
#endif

    typedef double Float;
    inline Float Lerp(Float t, Float v1, Float v2) { return (1 - t) * v1 + t * v2; }
    inline void CHECK_NE(Float a, Float b) {}
    inline void DCHECK(bool b) {}
    const Float Infinity = std::numeric_limits<Float>::infinity();

    inline void CHECK_GT(Float a, Float b) {}
    inline void CHECK_LT(Float a, Float b) {}
    inline void CHECK(bool b) {}
    inline void CHECK_GE(int a, int b) {}

    inline void CHECK_EQ(char a,char b) {}




    static Float ShadowEpsilon = 0.0001f;
    static Float Pi = 3.14159265358979323846;
    static Float InvPi = 0.31830988618379067154;
    static Float Inv2Pi = 0.15915494309189533577;
    static Float Inv4Pi = 0.07957747154594766788;
    static Float PiOver2 = 1.57079632679489661923;
    static Float PiOver4 = 0.78539816339744830961;
    static Float Sqrt2 = 1.41421356237309504880;




    template <typename T, typename U, typename V>
    inline T Clamp(T val, U low, V high) {
        if (val < low)
            return low;
        else if (val > high)
            return high;
        else
            return val;
    }






}



#endif