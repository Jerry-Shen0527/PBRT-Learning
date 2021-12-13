#ifndef PDF_H
#define PDF_H
#include "vec3.h"
#include "onb.h"

class pdf {
public:
    virtual ~pdf() {}

    virtual double value(const Vector3f& direction) const = 0;
    virtual Vector3f generate() const = 0;
};

class cosine_pdf : public pdf {
public:
    cosine_pdf(const Vector3f& w) { uvw.build_from_w(w); }
    cosine_pdf(const Normal3f& w) { uvw.build_from_w(w); }
    //caculate sample_pdf
    virtual double value(const Vector3f& direction) const override {
        auto cosine = Dot(unit_vector(direction), uvw.w());
        return (cosine <= 0) ? 0 : cosine / pi;
    }

    //generate/sample a direction according some distribution like random_cosine_direction()
    virtual Vector3f generate() const override {
        return uvw.local(random_cosine_direction());
    }

public:
    onb uvw;
};


//sample in hittable like plane light
class hittable_pdf : public pdf {
public:
    hittable_pdf(shared_ptr<hittable> p, const Point3f& origin) : ptr(p), o(origin) {}

    //caculate sample_pdf
    virtual double value(const Vector3f& direction) const override {
        return ptr->pdf_value(o, direction);
    }

    //generate/sample a direction according some distribution like random_cosine_direction()
    virtual Vector3f generate() const override {
        return ptr->random(o);
    }

public:
    Point3f o;
    shared_ptr<hittable> ptr;
};

class mixture_pdf : public pdf {
public:
    mixture_pdf(shared_ptr<pdf> p0, shared_ptr<pdf> p1) {
        p[0] = p0;
        p[1] = p1;
    }

    virtual double value(const Vector3f& direction) const override {
        return 0.5 * p[0]->value(direction) + 0.5 * p[1]->value(direction);
    }

    virtual Vector3f generate() const override {
        if (random_Float() < 0.5)
            return p[0]->generate();
        else
            return p[1]->generate();
    }

public:
    shared_ptr<pdf> p[2];
};

inline Vector3f random_to_sphere(double radius, double distance_squared) {
    auto r1 = random_Float();
    auto r2 = random_Float();
    auto z = 1 + r2 * (sqrt(1 - radius * radius / distance_squared) - 1);

    auto phi = 2 * pi * r1;
    auto x = cos(phi) * sqrt(1 - z * z);
    auto y = sin(phi) * sqrt(1 - z * z);

    return Vector3f(x, y, z);
}

#endif