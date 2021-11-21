#pragma once
#ifndef SPECTRUM
#define SPECTRUM

#include <iostream>
#include "rtweekend.h"

#define CHECKNAN(condition);
#define DCHECK(condition);
#define CHECKSORTED(lambda,vals,n);

#ifdef DEBUG
#define CHECKNAN(condition) if(condition) std::cout << "exist NaN!" << std::endl;
#define DCHECK(condition) if(!condition) std::cout<<"out of range!"<<std::endl;
#define CHECKSORTED(lambda,vals,n) if(!SpectrumSamplesSorted(lambda,vals,n)) std::cout<<"is not sorted!"<<std::endl;
#endif // DEBUG


static const int sampledLambdaStart = 400;
static const int sampledLambdaEnd = 700;
static const int nSpectralSamples = 60;

extern bool SpectrumSamplesSorted(const Float * lambda, const Float * vals,
    int n);
extern void SortSpectrumSamples(Float * lambda, Float * vals, int n);
extern Float AverageSpectrumSamples(const Float * lambda, const Float * vals,
    int n, Float lambdaStart, Float lambdaEnd);
inline void XYZToRGB(const Float xyz[3], Float * rgb) {
    rgb[0] = 3.240479f * xyz[0] - 1.537150f * xyz[1] - 0.498535f * xyz[2];
    rgb[1] = -0.969256f * xyz[0] + 1.875991f * xyz[1] + 0.041556f * xyz[2];
    rgb[2] = 0.055648f * xyz[0] - 0.204043f * xyz[1] + 1.057311f * xyz[2];
}
inline void XYZToRGB(const Float xyz[3], color& rgb) {
    Float frgb[3];
    XYZToRGB(xyz, frgb);
    rgb[0] = frgb[0]; rgb[1] = frgb[1]; rgb[2] = frgb[2];
}

inline void RGBToXYZ(const Float rgb[3], Float xyz[3]) {
    xyz[0] = 0.412453f * rgb[0] + 0.357580f * rgb[1] + 0.180423f * rgb[2];
    xyz[1] = 0.212671f * rgb[0] + 0.715160f * rgb[1] + 0.072169f * rgb[2];
    xyz[2] = 0.019334f * rgb[0] + 0.119193f * rgb[1] + 0.950227f * rgb[2];
}

enum class SpectrumType { Reflectance, Illuminant };
extern Float InterpolateSpectrumSamples(const Float * lambda, const Float * vals,
    int n, Float l);
extern void Blackbody(const Float * lambda, int n, Float T, Float * Le);
extern void BlackbodyNormalized(const Float * lambda, int n, Float T,
    Float * vals);

// Spectral Data Declarations
static const int nCIESamples = 471;
extern const Float CIE_X[nCIESamples];
extern const Float CIE_Y[nCIESamples];
extern const Float CIE_Z[nCIESamples];
extern const Float CIE_lambda[nCIESamples];
static const Float CIE_Y_integral = 106.856895;
static const int nRGB2SpectSamples = 32;
extern const Float RGB2SpectLambda[nRGB2SpectSamples];
extern const Float RGBRefl2SpectWhite[nRGB2SpectSamples];
extern const Float RGBRefl2SpectCyan[nRGB2SpectSamples];
extern const Float RGBRefl2SpectMagenta[nRGB2SpectSamples];
extern const Float RGBRefl2SpectYellow[nRGB2SpectSamples];
extern const Float RGBRefl2SpectRed[nRGB2SpectSamples];
extern const Float RGBRefl2SpectGreen[nRGB2SpectSamples];
extern const Float RGBRefl2SpectBlue[nRGB2SpectSamples];
extern const Float RGBIllum2SpectWhite[nRGB2SpectSamples];
extern const Float RGBIllum2SpectCyan[nRGB2SpectSamples];
extern const Float RGBIllum2SpectMagenta[nRGB2SpectSamples];
extern const Float RGBIllum2SpectYellow[nRGB2SpectSamples];
extern const Float RGBIllum2SpectRed[nRGB2SpectSamples];
extern const Float RGBIllum2SpectGreen[nRGB2SpectSamples];
extern const Float RGBIllum2SpectBlue[nRGB2SpectSamples];


template <int nSpectrumSamples>
class CoefficientSpectrum {
public:
    // CoefficientSpectrum Public Methods
    CoefficientSpectrum(Float v = 0.f) {
        for (int i = 0; i < nSpectrumSamples; i++) c[i] = v;
        CHECKNAN(HasNaNs());
    }

    CoefficientSpectrum operator+(const CoefficientSpectrum& v2) const
    {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; i++)
            ret.c[i] += v2.c[i];
        CHECKNAN(ret.HasNaNs());
        return ret;
    }

    CoefficientSpectrum& operator+=(const CoefficientSpectrum& v2)
    {
        for (int i = 0; i < nSpectrumSamples; i++)
            c[i] += v2.c[i];
        CHECKNAN(HasNaNs());
        return *this;
    }

    CoefficientSpectrum operator-(const CoefficientSpectrum& v2) const
    {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; i++)
            ret.c[i] -= v2.c[i];
        CHECKNAN(ret.HasNaNs());
        return ret;
    }

    CoefficientSpectrum& operator-=(const CoefficientSpectrum& v2)
    {
        for (int i = 0; i < nSpectrumSamples; i++)
            c[i] -= v2.c[i];
        CHECKNAN(HasNaNs());
        return *this;
    }

    CoefficientSpectrum operator*(const CoefficientSpectrum& v2) const
    {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; i++)
            ret.c[i] *= v2.c[i];
        CHECKNAN(ret.HasNaNs());
        return ret;
    }

    CoefficientSpectrum& operator*=(const CoefficientSpectrum& v2)
    {
        for (int i = 0; i < nSpectrumSamples; i++)
            c[i] *= v2.c[i];
        CHECKNAN(HasNaNs());
        return *this;
    }

    CoefficientSpectrum operator/(const CoefficientSpectrum& v2) const
    {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; i++)
            ret.c[i] /= v2.c[i];
        CHECKNAN(ret.HasNaNs());
        return ret;
    }

    CoefficientSpectrum& operator/=(const CoefficientSpectrum& v2)
    {
        for (int i = 0; i < nSpectrumSamples; i++)
            c[i] /= v2.c[i];
        CHECKNAN(HasNaNs());
        return *this;
    }

    CoefficientSpectrum operator*(Float f) const
    {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; i++)
            ret.c[i] *= f;
        CHECKNAN(ret.HasNaNs());
        return ret;
    }

    CoefficientSpectrum& operator*=(Float f)
    {
        for (int i = 0; i < nSpectrumSamples; i++)
            c[i] *= f;
        CHECKNAN(HasNaNs());
        return *this;
    }

    friend inline CoefficientSpectrum operator*(Float f, const CoefficientSpectrum& v)
    {
        return v * f;
    }

    CoefficientSpectrum operator/(Float f) const
    {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; i++)
            ret.c[i] /= f;
        CHECKNAN(ret.HasNaNs());
        return ret;
    }

    CoefficientSpectrum &operator/=(Float f) 
    {
        for (int i = 0; i < nSpectrumSamples; i++)
            c[i] /= f;
        CHECKNAN(HasNaNs());
        return ret;
    }

    CoefficientSpectrum operator-() const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] = -c[i];
        return ret;
    }

    bool operator==(const CoefficientSpectrum& sp) const {
        for (int i = 0; i < nSpectrumSamples; ++i)
            if (c[i] != sp.c[i]) return false;
        return true;
    }

    bool operator!=(const CoefficientSpectrum& sp) const {
        return !(*this == sp);
    }

    bool IsBlack() const {
        for (int i = 0; i < nSpectrumSamples; ++i)
            if (c[i] != 0.) return false;
        return true;
    }

    friend CoefficientSpectrum Sqrt(const CoefficientSpectrum& s) {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] = std::sqrt(s.c[i]);
        CHECKNAN(ret.HasNaNs());
        return ret;
    }

    friend std::ostream& operator<<(std::ostream& os, const CoefficientSpectrum& s) {
        os << "[";
        for (int i = 0; i < nSpectrumSamples; i++)
            os << s.c[i] << "\n";
        os << "]";
        return os;
    }

    CoefficientSpectrum Clamp(Float low = 0, Float high = infinity) const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSpectrumSamples; ++i)
            ret.c[i] = clamp(c[i], low, high);
        CHECKNAN(ret.HasNaNs());
        return ret;
    }

    Float MaxComponentValue() const {
        Float m = c[0];
        for (int i = 1; i < nSpectrumSamples; ++i)
            m = std::max(m, c[i]);
        return m;
    }

    bool Write(FILE* f) const {
        for (int i = 0; i < nSpectrumSamples; ++i)
            if (fprintf(f, "%f ", c[i]) < 0) return false;
        return true;
    }
    bool Read(FILE* f) {
        for (int i = 0; i < nSpectrumSamples; ++i) {
            double v;
            if (fscanf(f, "%lf ", &v) != 1) return false;
            c[i] = v;
        }
        return true;
    }
    Float& operator[](int i) {
        DCHECK(i >= 0 && i < nSpectrumSamples);
        return c[i];
    }
    Float operator[](int i) const {
        DCHECK(i >= 0 && i < nSpectrumSamples);
        return c[i];
    }


    template <int n> friend inline CoefficientSpectrum<n> Pow(const CoefficientSpectrum<n>& s,
        Float e);
    template <int n> friend inline CoefficientSpectrum<n> Exp(const CoefficientSpectrum<n>& s);

    bool HasNaNs() const {
     for (int i = 0; i < nSpectrumSamples; ++i)
            if (std::isnan(c[i])) return true;
        return false;
    }

    // CoefficientSpectrum Public Data
    static const int nSamples = nSpectrumSamples;

    
protected:
    // CoefficientSpectrum Protected Data
    Float c[nSpectrumSamples];
};


class RGBSpectrum : public CoefficientSpectrum<3> {
    using CoefficientSpectrum<3>::c;

public:
    // RGBSpectrum Public Methods
    RGBSpectrum(Float v = 0.f) : CoefficientSpectrum<3>(v) {}
    RGBSpectrum(const CoefficientSpectrum<3>& v) : CoefficientSpectrum<3>(v) {}
    RGBSpectrum(const RGBSpectrum& s,
        SpectrumType type = SpectrumType::Reflectance) {
        *this = s;
    }
    static RGBSpectrum FromRGB(const Float rgb[3],
        SpectrumType type = SpectrumType::Reflectance) {
        RGBSpectrum s;
        s.c[0] = rgb[0];
        s.c[1] = rgb[1];
        s.c[2] = rgb[2];
        CHECKNAN(s.HasNaNs());
        return s;
    }
    static RGBSpectrum FromRGB(const color cc,
        SpectrumType type = SpectrumType::Reflectance) {
        Float rgb[3] = { cc[0],cc[1],cc[2] };
        return FromRGB(rgb, type);
    }

    void ToRGB(Float* rgb) const {
        rgb[0] = c[0];
        rgb[1] = c[1];
        rgb[2] = c[2];
    }
    void ToRGB(color& cc) const {
        cc[0] = c[0]; cc[1] = c[1]; cc[2] = c[2];
    }

    const RGBSpectrum& ToRGBSpectrum() const { return *this; }
    void ToXYZ(Float xyz[3]) const { RGBToXYZ(c, xyz); }
    static RGBSpectrum FromXYZ(const Float xyz[3],
        SpectrumType type = SpectrumType::Reflectance) {
        RGBSpectrum r;
        XYZToRGB(xyz, r.c);
        return r;
    }
    Float y() const {
        const Float YWeight[3] = { 0.212671f, 0.715160f, 0.072169f };
        return YWeight[0] * c[0] + YWeight[1] * c[1] + YWeight[2] * c[2];
    }
    static RGBSpectrum FromSampled(const Float* lambda, const Float* v, int n) {
        // Sort samples if unordered, use sorted for returned spectrum
        if (!SpectrumSamplesSorted(lambda, v, n)) {
            std::vector<Float> slambda(&lambda[0], &lambda[n]);
            std::vector<Float> sv(&v[0], &v[n]);
            SortSpectrumSamples(&slambda[0], &sv[0], n);
            return FromSampled(&slambda[0], &sv[0], n);
        }
        Float xyz[3] = { 0, 0, 0 };
        for (int i = 0; i < nCIESamples; ++i) {
            Float val = InterpolateSpectrumSamples(lambda, v, n, CIE_lambda[i]);
            xyz[0] += val * CIE_X[i];
            xyz[1] += val * CIE_Y[i];
            xyz[2] += val * CIE_Z[i];
        }
        Float scale = Float(CIE_lambda[nCIESamples - 1] - CIE_lambda[0]) /
            Float(CIE_Y_integral * nCIESamples);
        xyz[0] *= scale;
        xyz[1] *= scale;
        xyz[2] *= scale;
        return FromXYZ(xyz);
    }
};

class SampledSpectrum : public CoefficientSpectrum<nSpectralSamples>{
  public:
      // SampledSpectrum Public Methods
      SampledSpectrum(Float v = 0.f) : CoefficientSpectrum(v) {}
      SampledSpectrum(const CoefficientSpectrum<nSpectralSamples> &v)
          : CoefficientSpectrum<nSpectralSamples>(v) {}
      static SampledSpectrum FromSampled(const Float * lambda, const Float * v,
                                         int n) {
          // Sort samples if unordered, use sorted for returned spectrum
          if (!SpectrumSamplesSorted(lambda, v, n)) {
              std::vector<Float> slambda(&lambda[0], &lambda[n]);
              std::vector<Float> sv(&v[0], &v[n]);
              SortSpectrumSamples(&slambda[0], &sv[0], n);
              return FromSampled(&slambda[0], &sv[0], n);
          }
          SampledSpectrum r;
          for (int i = 0; i < nSpectralSamples; ++i) {
              // Compute average value of given SPD over $i$th sample's range
              Float lambda0 = Lerp(Float(i) / Float(nSpectralSamples),
                                   sampledLambdaStart, sampledLambdaEnd);
              Float lambda1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
                                   sampledLambdaStart, sampledLambdaEnd);
              r.c[i] = AverageSpectrumSamples(lambda, v, n, lambda0, lambda1);
          }
          return r;
      }

      static void Init() {
          // Compute XYZ matching functions for _SampledSpectrum_
          for (int i = 0; i < nSpectralSamples; ++i) {
              Float wl0 = Lerp(Float(i) / Float(nSpectralSamples),
                               sampledLambdaStart, sampledLambdaEnd);
              Float wl1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
                               sampledLambdaStart, sampledLambdaEnd);
              X.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_X, nCIESamples, wl0,
                                              wl1);
              Y.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Y, nCIESamples, wl0,
                                              wl1);
              Z.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Z, nCIESamples, wl0,
                                              wl1);
          }

          // Compute RGB to spectrum functions for _SampledSpectrum_
          for (int i = 0; i < nSpectralSamples; ++i) {
              Float wl0 = Lerp(Float(i) / Float(nSpectralSamples),
                               sampledLambdaStart, sampledLambdaEnd);
              Float wl1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
                               sampledLambdaStart, sampledLambdaEnd);
              rgbRefl2SpectWhite.c[i] =
                  AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectWhite,
                                         nRGB2SpectSamples, wl0, wl1);
              rgbRefl2SpectCyan.c[i] =
                  AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectCyan,
                                         nRGB2SpectSamples, wl0, wl1);
              rgbRefl2SpectMagenta.c[i] =
                  AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectMagenta,
                                         nRGB2SpectSamples, wl0, wl1);
              rgbRefl2SpectYellow.c[i] =
                  AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectYellow,
                                         nRGB2SpectSamples, wl0, wl1);
              rgbRefl2SpectRed.c[i] = AverageSpectrumSamples(
                  RGB2SpectLambda, RGBRefl2SpectRed, nRGB2SpectSamples, wl0, wl1);
              rgbRefl2SpectGreen.c[i] =
                  AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectGreen,
                                         nRGB2SpectSamples, wl0, wl1);
              rgbRefl2SpectBlue.c[i] =
                  AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectBlue,
                                         nRGB2SpectSamples, wl0, wl1);

              rgbIllum2SpectWhite.c[i] =
                  AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectWhite,
                                         nRGB2SpectSamples, wl0, wl1);
              rgbIllum2SpectCyan.c[i] =
                  AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectCyan,
                                         nRGB2SpectSamples, wl0, wl1);
              rgbIllum2SpectMagenta.c[i] =
                  AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectMagenta,
                                         nRGB2SpectSamples, wl0, wl1);
              rgbIllum2SpectYellow.c[i] =
                  AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectYellow,
                                         nRGB2SpectSamples, wl0, wl1);
              rgbIllum2SpectRed.c[i] =
                  AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectRed,
                                         nRGB2SpectSamples, wl0, wl1);
              rgbIllum2SpectGreen.c[i] =
                  AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectGreen,
                                         nRGB2SpectSamples, wl0, wl1);
              rgbIllum2SpectBlue.c[i] =
                  AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectBlue,
                                         nRGB2SpectSamples, wl0, wl1);
          }
      }

      void ToXYZ(Float xyz[3]) const {
          xyz[0] = xyz[1] = xyz[2] = 0.f;
          for (int i = 0; i < nSpectralSamples; ++i) {
              xyz[0] += X.c[i] * c[i];
              xyz[1] += Y.c[i] * c[i];
              xyz[2] += Z.c[i] * c[i];
          }
          Float scale = Float(sampledLambdaEnd - sampledLambdaStart) /
                        Float(CIE_Y_integral * nSpectralSamples);
          xyz[0] *= scale;
          xyz[1] *= scale;
          xyz[2] *= scale;
      }

      Float y() const {
          Float yy = 0.f;
          for (int i = 0; i < nSpectralSamples; ++i) yy += Y.c[i] * c[i];
          return yy * Float(sampledLambdaEnd - sampledLambdaStart) /
                 Float(CIE_Y_integral * nSpectralSamples);
      }
 
      void ToRGB(Float* rgb) const {
          Float xyz[3];
          ToXYZ(xyz);
          XYZToRGB(xyz, rgb);
      }
      void ToRGB(color& cc) const {
          Float xyz[3];
          ToXYZ(xyz);
          XYZToRGB(xyz, cc);
      }

      RGBSpectrum ToRGBSpectrum() const;
      static SampledSpectrum FromRGB(
          const Float rgb[3], SpectrumType type = SpectrumType::Illuminant);
      static SampledSpectrum FromRGB(
          const color rgb, SpectrumType type = SpectrumType::Illuminant) {
          Float frgb[3] = { rgb[0],rgb[1],rgb[2] };
          return FromRGB(frgb, type);
      }

      static SampledSpectrum FromXYZ(
          const Float xyz[3], SpectrumType type = SpectrumType::Reflectance) {
          Float rgb[3];
          XYZToRGB(xyz, rgb);
          return FromRGB(rgb, type);
      }

      SampledSpectrum(const RGBSpectrum& r,
                      SpectrumType type = SpectrumType::Reflectance);

    private:
        // SampledSpectrum Private Data
        static SampledSpectrum X, Y, Z;
        static SampledSpectrum rgbRefl2SpectWhite, rgbRefl2SpectCyan;
        static SampledSpectrum rgbRefl2SpectMagenta, rgbRefl2SpectYellow;
        static SampledSpectrum rgbRefl2SpectRed, rgbRefl2SpectGreen;
        static SampledSpectrum rgbRefl2SpectBlue;
        static SampledSpectrum rgbIllum2SpectWhite, rgbIllum2SpectCyan;
        static SampledSpectrum rgbIllum2SpectMagenta, rgbIllum2SpectYellow;
        static SampledSpectrum rgbIllum2SpectRed, rgbIllum2SpectGreen;
        static SampledSpectrum rgbIllum2SpectBlue;
};


template <int nSpectrumSamples>
inline CoefficientSpectrum<nSpectrumSamples> Pow(
    const CoefficientSpectrum<nSpectrumSamples>& s, Float e) {
    CoefficientSpectrum<nSpectrumSamples> ret;
    for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] = std::pow(s.c[i], e);
    DCHECK(!ret.HasNaNs());
    return ret;
}

template <int n>
inline CoefficientSpectrum<n> Exp(const CoefficientSpectrum<n>& s) {
    CoefficientSpectrum<n> ret;
    for (int i = 0; i < n; ++i) ret.c[i] = std::exp(s.c[i]);
    CHECKNAN(ret.HasNaNs());
    std::cout << "outside" << std::endl;
    return ret;
}

#endif // !SPECTRUM