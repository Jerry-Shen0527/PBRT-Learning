#pragma once
#include <iostream>
#include <vector>

#include "check.h"
#include "vec3.h"


static const int sampledLambdaStart = 400;
static const int sampledLambdaEnd = 700;
static const int nSpectralSamples = 60;
extern bool SpectrumSamplesSorted(const double* lambda, const double* vals,
    int n);
extern void SortSpectrumSamples(double* lambda, double* vals, int n);
extern double AverageSpectrumSamples(const double* lambda, const double* vals,
    int n, double lambdaStart, double lambdaEnd);
inline void XYZToRGB(const double xyz[3], double rgb[3]) {
    rgb[0] = 3.240479f * xyz[0] - 1.537150f * xyz[1] - 0.498535f * xyz[2];
    rgb[1] = -0.969256f * xyz[0] + 1.875991f * xyz[1] + 0.041556f * xyz[2];
    rgb[2] = 0.055648f * xyz[0] - 0.204043f * xyz[1] + 1.057311f * xyz[2];
}

inline void RGBToXYZ(const double rgb[3], double xyz[3]) {
    xyz[0] = 0.412453f * rgb[0] + 0.357580f * rgb[1] + 0.180423f * rgb[2];
    xyz[1] = 0.212671f * rgb[0] + 0.715160f * rgb[1] + 0.072169f * rgb[2];
    xyz[2] = 0.019334f * rgb[0] + 0.119193f * rgb[1] + 0.950227f * rgb[2];
}



enum class SpectrumType { Reflectance, Illuminant };
extern double InterpolateSpectrumSamples(const double* lambda, const double* vals,
    int n, double l);
extern void Blackbody(const double* lambda, int n, double T, double* Le);
extern void BlackbodyNormalized(const double* lambda, int n, double T,
    double* vals);

// Spectral Data Declarations
static const int nCIESamples = 471;
extern const double CIE_X[nCIESamples];
extern const double CIE_Y[nCIESamples];
extern const double CIE_Z[nCIESamples];
extern const double CIE_lambda[nCIESamples];
static const double CIE_Y_integral = 106.856895;
static const int nRGB2SpectSamples = 32;
extern const double RGB2SpectLambda[nRGB2SpectSamples];
extern const double RGBRefl2SpectWhite[nRGB2SpectSamples];
extern const double RGBRefl2SpectCyan[nRGB2SpectSamples];
extern const double RGBRefl2SpectMagenta[nRGB2SpectSamples];
extern const double RGBRefl2SpectYellow[nRGB2SpectSamples];
extern const double RGBRefl2SpectRed[nRGB2SpectSamples];
extern const double RGBRefl2SpectGreen[nRGB2SpectSamples];
extern const double RGBRefl2SpectBlue[nRGB2SpectSamples];
extern const double RGBIllum2SpectWhite[nRGB2SpectSamples];
extern const double RGBIllum2SpectCyan[nRGB2SpectSamples];
extern const double RGBIllum2SpectMagenta[nRGB2SpectSamples];
extern const double RGBIllum2SpectYellow[nRGB2SpectSamples];
extern const double RGBIllum2SpectRed[nRGB2SpectSamples];
extern const double RGBIllum2SpectGreen[nRGB2SpectSamples];
extern const double RGBIllum2SpectBlue[nRGB2SpectSamples];

template <int nSpectrumSamples>
class CoefficientSpectrum {
public:
    // CoefficientSpectrum Public Methods
    CoefficientSpectrum(double v = 0.f) {
        for (int i = 0; i < nSpectrumSamples; ++i) c[i] = v;
        check::DCHECK(!HasNaNs());
    }
#ifdef DEBUG
    CoefficientSpectrum(const CoefficientSpectrum& s) {
        DCHECK(!s.HasNaNs());
        for (int i = 0; i < nSpectrumSamples; ++i) c[i] = s.c[i];
    }

    CoefficientSpectrum& operator=(const CoefficientSpectrum& s) {
        DCHECK(!s.HasNaNs());
        for (int i = 0; i < nSpectrumSamples; ++i) c[i] = s.c[i];
        return *this;
    }
#endif  // DEBUG
    void Print(FILE* f) const {
        fprintf(f, "[ ");
        for (int i = 0; i < nSpectrumSamples; ++i) {
            fprintf(f, "%f", c[i]);
            if (i != nSpectrumSamples - 1) fprintf(f, ", ");
        }
        fprintf(f, "]");
    }
    CoefficientSpectrum& operator+=(const CoefficientSpectrum& s2) {
        check::DCHECK(!s2.HasNaNs());
        for (int i = 0; i < nSpectrumSamples; ++i) c[i] += s2.c[i];
        return *this;
    }
    CoefficientSpectrum operator+(const CoefficientSpectrum& s2) const {
        check::DCHECK(!s2.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] += s2.c[i];
        return ret;
    }
    CoefficientSpectrum operator-(const CoefficientSpectrum& s2) const {
        check::DCHECK(!s2.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] -= s2.c[i];
        return ret;
    }
    CoefficientSpectrum operator/(const CoefficientSpectrum& s2) const {
        check::DCHECK(!s2.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) {
            CHECK_NE(s2.c[i], 0);
            ret.c[i] /= s2.c[i];
        }
        return ret;
    }
    CoefficientSpectrum operator*(const CoefficientSpectrum& sp) const {
        check::DCHECK(!sp.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] *= sp.c[i];
        return ret;
    }
    CoefficientSpectrum& operator*=(const CoefficientSpectrum& sp) {
        check::DCHECK(!sp.HasNaNs());
        for (int i = 0; i < nSpectrumSamples; ++i) c[i] *= sp.c[i];
        return *this;
    }
    CoefficientSpectrum operator*(double a) const {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] *= a;
        check::DCHECK(!ret.HasNaNs());
        return ret;
    }
    CoefficientSpectrum& operator*=(double a) {
        for (int i = 0; i < nSpectrumSamples; ++i) c[i] *= a;
        check::DCHECK(!HasNaNs());
        return *this;
    }
    friend inline CoefficientSpectrum operator*(double a,
        const CoefficientSpectrum& s) {
        check::DCHECK(!std::isnan(a) && !s.HasNaNs());
        return s * a;
    }
    CoefficientSpectrum operator/(double a) const {
        //check::CHECK_NE(a, 0);
        check::DCHECK(!std::isnan(a));
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] /= a;
        DCHECK(!ret.HasNaNs());
        return ret;
    }
    CoefficientSpectrum& operator/=(double a) {
        //CHECK_NE(a, 0);
        check::DCHECK(!std::isnan(a));
        for (int i = 0; i < nSpectrumSamples; ++i) c[i] /= a;
        return *this;
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
        check::DCHECK(!ret.HasNaNs());
        return ret;
    }
    template <int n>
    friend inline CoefficientSpectrum<n> Pow(const CoefficientSpectrum<n>& s,
        double e);
    CoefficientSpectrum operator-() const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] = -c[i];
        return ret;
    }
    friend CoefficientSpectrum Exp(const CoefficientSpectrum& s) {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] = std::exp(s.c[i]);
        check::DCHECK(!ret.HasNaNs());
        return ret;
    }
    friend std::ostream& operator<<(std::ostream& os,
        const CoefficientSpectrum& s) {
        return os << s.ToString();
    }
    std::string ToString() const {
        std::string str = "[ ";
        for (int i = 0; i < nSpectrumSamples; ++i) {
            str += StringPrintf("%f", c[i]);
            if (i + 1 < nSpectrumSamples) str += ", ";
        }
        str += " ]";
        return str;
    }
    CoefficientSpectrum Clamp(double low = 0, double high = INFINITY) const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSpectrumSamples; ++i)
            ret.c[i] = check::Clamp(c[i], low, high);
        check::DCHECK(!ret.HasNaNs());
        return ret;
    }
    double MaxComponentValue() const {
        double m = c[0];
        for (int i = 1; i < nSpectrumSamples; ++i)
            m = std::max(m, c[i]);
        return m;
    }
    bool HasNaNs() const {
        for (int i = 0; i < nSpectrumSamples; ++i)
            if (std::isnan(c[i])) return true;
        return false;
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
    double& operator[](int i) {
        check::DCHECK(i >= 0 && i < nSpectrumSamples);
        return c[i];
    }
    double operator[](int i) const {
        check::DCHECK(i >= 0 && i < nSpectrumSamples);
        return c[i];
    }

    // CoefficientSpectrum Public Data
    static const int nSamples = nSpectrumSamples;

protected:
    // CoefficientSpectrum Protected Data
    double c[nSpectrumSamples];
};

class SampledSpectrum : public CoefficientSpectrum<nSpectralSamples> {
public:
    // SampledSpectrum Public Methods
    SampledSpectrum(double v = 0.f) : CoefficientSpectrum(v) {}
    SampledSpectrum(const CoefficientSpectrum<nSpectralSamples>& v)
        : CoefficientSpectrum<nSpectralSamples>(v) {}
    /*SampledSpectrum(const RGBSpectrum& r,
        SpectrumType type = SpectrumType::Reflectance);*/


    static SampledSpectrum FromSampled(const double* lambda, const double* v,
        int n) {
        // Sort samples if unordered, use sorted for returned spectrum
        if (!SpectrumSamplesSorted(lambda, v, n)) {
            std::vector<double> slambda(&lambda[0], &lambda[n]);
            std::vector<double> sv(&v[0], &v[n]);
            SortSpectrumSamples(&slambda[0], &sv[0], n);
            return FromSampled(&slambda[0], &sv[0], n);
        }
        SampledSpectrum r;
        for (int i = 0; i < nSpectralSamples; ++i) {
            // Compute average value of given SPD over $i$th sample's range
            double lambda0 = check::Lerp(double(i) / double(nSpectralSamples),
                sampledLambdaStart, sampledLambdaEnd);
            double lambda1 = check::Lerp(double(i + 1) / double(nSpectralSamples),
                sampledLambdaStart, sampledLambdaEnd);
            r.c[i] = AverageSpectrumSamples(lambda, v, n, lambda0, lambda1);
        }
        return r;
    }
    static void Init() {
        // Compute XYZ matching functions for _SampledSpectrum_
        for (int i = 0; i < nSpectralSamples; ++i) {
            double wl0 = check::Lerp(double(i) / double(nSpectralSamples),
                sampledLambdaStart, sampledLambdaEnd);
            double wl1 = check::Lerp(double(i + 1) / double(nSpectralSamples),
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
            double wl0 = check::Lerp(double(i) / double(nSpectralSamples),
                sampledLambdaStart, sampledLambdaEnd);
            double wl1 = check::Lerp(double(i + 1) / double(nSpectralSamples),
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
    void ToXYZ(double xyz[3]) const {
        xyz[0] = xyz[1] = xyz[2] = 0.f;
        for (int i = 0; i < nSpectralSamples; ++i) {
            xyz[0] += X.c[i] * c[i];
            xyz[1] += Y.c[i] * c[i];
            xyz[2] += Z.c[i] * c[i];
        }
        double scale = double(sampledLambdaEnd - sampledLambdaStart) /
            double(CIE_Y_integral * nSpectralSamples);
        xyz[0] *= scale;
        xyz[1] *= scale;
        xyz[2] *= scale;
    }
    double y() const {
        double yy = 0.f;
        for (int i = 0; i < nSpectralSamples; ++i) yy += Y.c[i] * c[i];
        return yy * double(sampledLambdaEnd - sampledLambdaStart) /
            double(CIE_Y_integral * nSpectralSamples);
    }
    void ToRGB(double rgb[3]) const {
        double xyz[3];
        ToXYZ(xyz);
        XYZToRGB(xyz, rgb);
    }
    // RGBSpectrum ToRGBSpectrum() const;
    static SampledSpectrum FromRGB(
        const double rgb[3], SpectrumType type = SpectrumType::Illuminant);
    static SampledSpectrum FromXYZ(
        const double xyz[3], SpectrumType type = SpectrumType::Reflectance) {
        double rgb[3];
        XYZToRGB(xyz, rgb);
        return FromRGB(rgb, type);
    }


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

class RGBSpectrum : public CoefficientSpectrum<3> {
    using CoefficientSpectrum<3>::c;

public:
    // RGBSpectrum Public Methods
    RGBSpectrum(double v = 0.f) : CoefficientSpectrum<3>(v) {}
    RGBSpectrum(const CoefficientSpectrum<3>& v) : CoefficientSpectrum<3>(v) {}
    RGBSpectrum(const RGBSpectrum& s,
        SpectrumType type = SpectrumType::Reflectance) {
        *this = s;
    }
    static RGBSpectrum FromRGB(const double rgb[3],
        SpectrumType type = SpectrumType::Reflectance) {
        RGBSpectrum s;
        s.c[0] = rgb[0];
        s.c[1] = rgb[1];
        s.c[2] = rgb[2];
        //DCHECK(!s.HasNaNs());
        return s;
    }
    void ToRGB(double* rgb) const {
        rgb[0] = c[0];
        rgb[1] = c[1];
        rgb[2] = c[2];
    }
    color ToColor() const {
        color color_tmp;
        ToRGB(color_tmp.e);
        return color_tmp;
    }
    const RGBSpectrum& ToRGBSpectrum() const { return *this; }
    void ToXYZ(double xyz[3]) const { RGBToXYZ(c, xyz); }
    static RGBSpectrum FromXYZ(const double xyz[3],
        SpectrumType type = SpectrumType::Reflectance) {
        RGBSpectrum r;
        XYZToRGB(xyz, r.c);
        return r;
    }
    double y() const {
        const double YWeight[3] = { 0.212671f, 0.715160f, 0.072169f };
        return YWeight[0] * c[0] + YWeight[1] * c[1] + YWeight[2] * c[2];
    }
    static RGBSpectrum FromSampled(const double* lambda, const double* v, int n) {
        // Sort samples if unordered, use sorted for returned spectrum
        if (!SpectrumSamplesSorted(lambda, v, n)) {
            std::vector<double> slambda(&lambda[0], &lambda[n]);
            std::vector<double> sv(&v[0], &v[n]);
            SortSpectrumSamples(&slambda[0], &sv[0], n);
            return FromSampled(&slambda[0], &sv[0], n);
        }
        double xyz[3] = { 0, 0, 0 };
        for (int i = 0; i < nCIESamples; ++i) {
            double val = InterpolateSpectrumSamples(lambda, v, n, CIE_lambda[i]);
            xyz[0] += val * CIE_X[i];
            xyz[1] += val * CIE_Y[i];
            xyz[2] += val * CIE_Z[i];
        }
        double scale = double(CIE_lambda[nCIESamples - 1] - CIE_lambda[0]) /
            double(CIE_Y_integral * nCIESamples);
        xyz[0] *= scale;
        xyz[1] *= scale;
        xyz[2] *= scale;
        return FromXYZ(xyz);
    }
};

// Spectrum Inline Functions


void ResampleLinearSpectrum(const double* lambdaIn, const double* vIn, int nIn,
    double lambdaMin, double lambdaMax, int nOut,
    double* vOut);


