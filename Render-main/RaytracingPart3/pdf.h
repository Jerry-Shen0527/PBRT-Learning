#ifndef PDF_H
#define PDF_H
#include "geo.h"
#include "onb.h"
#include "Primitive.h"

class pdf {
public:
    virtual ~pdf() {}

    virtual double value(const Vector3f& direction) const = 0;
    virtual Vector3f generate() const = 0;
};

class cosine_pdf : public pdf {
public:
    cosine_pdf(const Vector3f& w) { uvw.build_from_w(w); }

    virtual double value(const Vector3f& direction) const override {
        auto cosine = Dot(Normalize(direction), uvw.w());
        //cout << direction.x << " " << direction.y << " " << endl;
        return (cosine <= 0) ? 0 : cosine / PI;
    }

    virtual Vector3f generate() const override {
        return uvw.local(random_cosine_direction());
    }

public:
    onb uvw;
};

class Primitive_pdf : public pdf {
public:
    Primitive_pdf(shared_ptr<Primitive> p, const Point3f& origin) : ptr(p), o(origin) {}

    virtual double value(const Vector3f& direction) const override {
        //cout << ptr->pdf_value(o, direction) << "Primitive" << endl;
        return ptr->pdf_value(o, direction);
    }

    virtual Vector3f generate() const override {
        return ptr->random(o);
    }

public:
    Point3f o;
    shared_ptr<Primitive> ptr;
};

class mixture_pdf : public pdf {
public:
    mixture_pdf(shared_ptr<pdf> p0, shared_ptr<pdf> p1) {
        p[0] = p0;
        p[1] = p1;
    }

    virtual double value(const Vector3f& direction) const override {
        //cout << p[0]->value(direction) << " " << p[1]->value(direction) << endl;
        //if(isNaN(p[0]->value(direction))) return p[1]->value(direction);
        return 0.5 * p[0]->value(direction) + 0.5 * p[1]->value(direction);
        //return p[1]->value(direction);
    }

    virtual Vector3f generate() const override {
        if (random_double() < 0.5)
            return p[0]->generate();
        else
            return p[1]->generate();
    }

public:
    shared_ptr<pdf> p[2];
};

#endif