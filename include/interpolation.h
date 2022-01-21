
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_INTERPOLATION_H
#define PBRT_CORE_INTERPOLATION_H

// core/interpolation.h*
//#include "pbrt.h"
#include "rtweekend.h"

namespace pbrt {

// Spline Interpolation Declarations
Float CatmullRom(int size, const Float *nodes, const Float *values, Float x);
bool CatmullRomWeights(int size, const Float *nodes, Float x, int *offset,
                       Float *weights);
Float SampleCatmullRom(int size, const Float *nodes, const Float *f,
                       const Float *cdf, Float sample, Float *fval = nullptr,
                       Float *pdf = nullptr);
Float SampleCatmullRom2D(int size1, int size2, const Float *nodes1,
                         const Float *nodes2, const Float *values,
                         const Float *cdf, Float alpha, Float sample,
                         Float *fval = nullptr, Float *pdf = nullptr);
Float IntegrateCatmullRom(int n, const Float *nodes, const Float *values,
                          Float *cdf);
Float InvertCatmullRom(int n, const Float *x, const Float *values, Float u);

// Fourier Interpolation Declarations
Float Fourier(const Float *a, int m, double cosPhi);
Float SampleFourier(const Float *ak, const Float *recip, int m, Float u,
                    Float *pdf, Float *phiPtr);

}  // namespace pbrt

#endif  // PBRT_CORE_INTERPOLATION_H
