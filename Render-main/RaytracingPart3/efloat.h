#pragma once
#include <iostream>
#include "check.h"
typedef signed char        int8_t;
typedef short              int16_t;
typedef int                int32_t;
typedef long long          int64_t;
typedef unsigned char      uint8_t;
typedef unsigned short     uint16_t;
typedef unsigned int       uint32_t;
typedef unsigned long long uint64_t;

using namespace std;
#ifdef PBRT_FLOAT_AS_DOUBLE
typedef double Float;
#else
typedef float Float;
#endif
#ifdef _MSC_VER
#define MachineEpsilon (std::numeric_limits<Float>::epsilon() * 0.5)
#else
static PBRT_CONSTEXPR Float MachineEpsilon =
std::numeric_limits<Float>::epsilon() * 0.5;
#endif


class EFloat {
public:
    // EFloat Public Methods
    EFloat() {}
    EFloat(float v, float err = 0.f) : v(v) {
        if (err == 0.)
            low = high = v;
        else {
            // Compute conservative bounds by rounding the endpoints away
            // from the middle. Note that this will be over-conservative in
            // cases where v-err or v+err are exactly representable in
            // floating-point, but it's probably not worth the trouble of
            // checking this case.
            low = NextFloatDown(v - err);
            high = NextFloatUp(v + err);
        }
        // Store high precision reference value in _EFloat_

        vPrecise = v;
        Check();

    }

    EFloat(float v, long double lD, float err) : EFloat(v, err) {
        vPrecise = lD;
        Check();
    }

    EFloat operator+(EFloat ef) const {
        EFloat r;
        r.v = v + ef.v;

        r.vPrecise = vPrecise + ef.vPrecise;
 // DEBUG
        // Interval arithemetic addition, with the result rounded away from
        // the value r.v in order to be conservative.
        r.low = NextFloatDown(LowerBound() + ef.LowerBound());
        r.high = NextFloatUp(UpperBound() + ef.UpperBound());
        r.Check();
        return r;
    }
    explicit operator float() const { return v; }
    explicit operator double() const { return v; }
    float GetAbsoluteError() const {
        return NextFloatUp(std::max(std::abs(high - v),
            std::abs(v - low)));
    }
    float UpperBound() const { return high; }
    float LowerBound() const { return low; }

    float GetRelativeError() const {
        return std::abs((vPrecise - v) / vPrecise);
    }
    long double PreciseValue() const { return vPrecise; }

    EFloat operator-(EFloat ef) const {
        EFloat r;
        r.v = v - ef.v;

        r.vPrecise = vPrecise - ef.vPrecise;

        r.low = NextFloatDown(LowerBound() - ef.UpperBound());
        r.high = NextFloatUp(UpperBound() - ef.LowerBound());
        r.Check();
        return r;
    }
    EFloat operator*(EFloat ef) const {
        EFloat r;
        r.v = v * ef.v;

        r.vPrecise = vPrecise * ef.vPrecise;

        Float prod[4] = {
            LowerBound() * ef.LowerBound(), UpperBound() * ef.LowerBound(),
            LowerBound() * ef.UpperBound(), UpperBound() * ef.UpperBound() };
        r.low = NextFloatDown(
            std::min(std::min(prod[0], prod[1]), std::min(prod[2], prod[3])));
        r.high = NextFloatUp(
            std::max(std::max(prod[0], prod[1]), std::max(prod[2], prod[3])));
        r.Check();
        return r;
    }
    EFloat operator/(EFloat ef) const {
        EFloat r;
        r.v = v / ef.v;

        r.vPrecise = vPrecise / ef.vPrecise;

        if (ef.low < 0 && ef.high > 0) {
            // Bah. The interval we're dividing by straddles zero, so just
            // return an interval of everything.
            r.low = -INFINITY;
            r.high = INFINITY;
        }
        else {
            Float div[4] = {
                LowerBound() / ef.LowerBound(), UpperBound() / ef.LowerBound(),
                LowerBound() / ef.UpperBound(), UpperBound() / ef.UpperBound() };
            r.low = NextFloatDown(
                std::min(std::min(div[0], div[1]), std::min(div[2], div[3])));
            r.high = NextFloatUp(
                std::max(std::max(div[0], div[1]), std::max(div[2], div[3])));
        }
        r.Check();
        return r;
    }
    EFloat operator-() const {
        EFloat r;
        r.v = -v;

        r.vPrecise = -vPrecise;

        r.low = -high;
        r.high = -low;
        r.Check();
        return r;
    }
    inline bool operator==(EFloat fe) const { return v == fe.v; }
    inline void Check() const {
        if (!std::isinf(low) && !std::isnan(low) && !std::isinf(high) &&
            !std::isnan(high))
            check::CHECK_LE(low, high);

        if (!std::isinf(v) && !std::isnan(v)) {
            check::CHECK_LE(LowerBound(), vPrecise);
            check::CHECK_LE(vPrecise, UpperBound());
        }

    }
    EFloat(const EFloat& ef) {
        ef.Check();
        v = ef.v;
        low = ef.low;
        high = ef.high;

        vPrecise = ef.vPrecise;

    }
    EFloat& operator=(const EFloat& ef) {
        ef.Check();
        if (&ef != this) {
            v = ef.v;
            low = ef.low;
            high = ef.high;

            vPrecise = ef.vPrecise;

        }
        return *this;
    }

    /*friend std::ostream& operator<<(std::ostream& os, const EFloat& ef) {
        os << printf("v=%f (%a) - [%f, %f]",
            ef.v, ef.v, ef.low, ef.high);

        os << StringPrintf(", precise=%.30Lf", ef.vPrecise);

        return os;
    }*/

private:
    // EFloat Private Data
    float v, low, high;

    long double vPrecise;

    friend inline EFloat sqrt(EFloat fe);
    friend inline EFloat abs(EFloat fe);
    friend inline bool Quadratic(EFloat A, EFloat B, EFloat C, EFloat* t0,
        EFloat* t1);
};

// EFloat Inline Functions
inline EFloat operator*(float f, EFloat fe) { return EFloat(f) * fe; }

inline EFloat operator/(float f, EFloat fe) { return EFloat(f) / fe; }

inline EFloat operator+(float f, EFloat fe) { return EFloat(f) + fe; }

inline EFloat operator-(float f, EFloat fe) { return EFloat(f) - fe; }

inline EFloat sqrt(EFloat fe) {
    EFloat r;
    r.v = std::sqrt(fe.v);

    r.vPrecise = std::sqrt(fe.vPrecise);

    r.low = NextFloatDown(std::sqrt(fe.low));
    r.high = NextFloatUp(std::sqrt(fe.high));
    r.Check();
    return r;
}

inline EFloat abs(EFloat fe) {
    if (fe.low >= 0)
        // The entire interval is greater than zero, so we're all set.
        return fe;
    else if (fe.high <= 0) {
        // The entire interval is less than zero.
        EFloat r;
        r.v = -fe.v;

        r.vPrecise = -fe.vPrecise;

        r.low = -fe.high;
        r.high = -fe.low;
        r.Check();
        return r;
    }
    else {
        // The interval straddles zero.
        EFloat r;
        r.v = std::abs(fe.v);

        r.vPrecise = std::abs(fe.vPrecise);

        r.low = 0;
        r.high = std::max(-fe.low, fe.high);
        r.Check();
        return r;
    }
}

inline bool Quadratic(EFloat A, EFloat B, EFloat C, EFloat* t0, EFloat* t1);
inline bool Quadratic(EFloat A, EFloat B, EFloat C, EFloat* t0, EFloat* t1) {
    // Find quadratic discriminant
    double discrim = (double)B.v * (double)B.v - 4. * (double)A.v * (double)C.v;
    if (discrim < 0.) return false;
    double rootDiscrim = std::sqrt(discrim);

    EFloat floatRootDiscrim(rootDiscrim, MachineEpsilon * rootDiscrim);

    // Compute quadratic _t_ values
    EFloat q;
    if ((float)B < 0)
        q = -.5 * (B - floatRootDiscrim);
    else
        q = -.5 * (B + floatRootDiscrim);
    *t0 = q / A;
    *t1 = C / q;
    if ((float)*t0 > (float)*t1) std::swap(*t0, *t1);
    return true;
}