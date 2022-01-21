#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_TEXTURE_H
#define PBRT_CORE_TEXTURE_H

#include "rtweekend.h"
#include "color.h"
#include "perlin.h"
#include <iostream>

#include "transform.h"

class texture {
public:
    virtual color value(double u, double v, const Point3f& p) const = 0;
};

class solid_color : public texture {
public:
    solid_color() {}
    solid_color(color c) : color_value(c) {}

    solid_color(double red, double green, double blue)
        : solid_color(color(red, green, blue)) {}

    virtual color value(double u, double v, const Point3f& p) const override {
        return color_value;
    }

private:
    color color_value;
};

class checker_texture : public texture {
public:
    checker_texture() {}

    checker_texture(shared_ptr<texture> _even, shared_ptr<texture> _odd)
        : even(_even), odd(_odd) {}

    checker_texture(color c1, color c2)
        : even(make_shared<solid_color>(c1)), odd(make_shared<solid_color>(c2)) {}

    virtual color value(double u, double v, const Point3f& p) const override {
        auto sines = sin(10 * p.x) * sin(10 * p.y) * sin(10 * p.z);
        if (sines < 0)
            return odd->value(u, v, p);
        else
            return even->value(u, v, p);
    }

public:
    shared_ptr<texture> odd;
    shared_ptr<texture> even;
};

class noise_texture : public texture {
public:
    noise_texture() {}
    noise_texture(double sc) : scale(sc) {}

    virtual color value(double u, double v, const Point3f& p) const override {
        return color(1, 1, 1) * 0.5 * (1 + sin(scale * p.z + 10 * noise.turb(p)));
    }

public:
    perlin noise;
    double scale;
};

class image_texture : public texture {
public:
    const static int bytes_per_pixel = 3;

    image_texture()
        : data(nullptr), width(0), height(0), bytes_per_scanline(0) {}

    image_texture(const char* filename);

    ~image_texture() {
        delete data;
    }

    virtual color value(double u, double v, const Point3f& p) const override {
        // If we have no texture data, then return solid cyan as a debugging aid.
        
        if (data == nullptr)
            return color(0, 1, 1);
        
        // Clamp input texture coordinates to [0,1] x [1,0]
        u = clamp(u, 0.0, 1.0);
        v = 1.0 - clamp(v, 0.0, 1.0);  // Flip V to image coordinates

        auto i = static_cast<int>(u * width);
        auto j = static_cast<int>(v * height);

        // Clamp integer mapping, since actual coordinates should be less than 1.0
        if (i >= width)  i = width - 1;
        if (j >= height) j = height - 1;

        const auto color_scale = 1.0 / 255.0;
        auto pixel = data + j * bytes_per_scanline + i * bytes_per_pixel;

        return color(color_scale * pixel[0], color_scale * pixel[1], color_scale * pixel[2]);;
    }

private:
    unsigned char* data;
    int width, height;
    int bytes_per_scanline;
};

namespace pbrt {

    // Texture Declarations
    class TextureMapping2D {
    public:
        // TextureMapping2D Interface
        virtual ~TextureMapping2D();
        virtual Point2f Map(const SurfaceInteraction& si, Vector2f* dstdx,
            Vector2f* dstdy) const = 0;
    };

    class UVMapping2D : public TextureMapping2D {
    public:
        // UVMapping2D Public Methods
        UVMapping2D(Float su = 1, Float sv = 1, Float du = 0, Float dv = 0);
        Point2f Map(const SurfaceInteraction& si, Vector2f* dstdx,
            Vector2f* dstdy) const;

    private:
        const Float su, sv, du, dv;
    };

    class SphericalMapping2D : public TextureMapping2D {
    public:
        // SphericalMapping2D Public Methods
        SphericalMapping2D(const Transform& WorldToTexture)
            : WorldToTexture(WorldToTexture) {}
        Point2f Map(const SurfaceInteraction& si, Vector2f* dstdx,
            Vector2f* dstdy) const;

    private:
        Point2f sphere(const Point3f& P) const;
        const Transform WorldToTexture;
    };

    class CylindricalMapping2D : public TextureMapping2D {
    public:
        // CylindricalMapping2D Public Methods
        CylindricalMapping2D(const Transform& WorldToTexture)
            : WorldToTexture(WorldToTexture) {}
        Point2f Map(const SurfaceInteraction& si, Vector2f* dstdx,
            Vector2f* dstdy) const;

    private:
        // CylindricalMapping2D Private Methods
        Point2f cylinder(const Point3f& p) const {
            Vector3f vec = Normalize(WorldToTexture(p) - Point3f(0, 0, 0));
            return Point2f((Pi + std::atan2(vec.y, vec.x)) * Inv2Pi, vec.z);
        }
        const Transform WorldToTexture;
    };

    class PlanarMapping2D : public TextureMapping2D {
    public:
        // PlanarMapping2D Public Methods
        Point2f Map(const SurfaceInteraction& si, Vector2f* dstdx,
            Vector2f* dstdy) const;
        PlanarMapping2D(const Vector3f& vs, const Vector3f& vt, Float ds = 0,
            Float dt = 0)
            : vs(vs), vt(vt), ds(ds), dt(dt) {}

    private:
        const Vector3f vs, vt;
        const Float ds, dt;
    };

    class TextureMapping3D {
    public:
        // TextureMapping3D Interface
        virtual ~TextureMapping3D();
        virtual Point3f Map(const SurfaceInteraction& si, Vector3f* dpdx,
            Vector3f* dpdy) const = 0;
    };

    class IdentityMapping3D : public TextureMapping3D {
    public:
        // IdentityMapping3D Public Methods
        IdentityMapping3D(const Transform& WorldToTexture)
            : WorldToTexture(WorldToTexture) {}
        Point3f Map(const SurfaceInteraction& si, Vector3f* dpdx,
            Vector3f* dpdy) const;

    private:
        const Transform WorldToTexture;
    };

    template <typename T>
    class Texture {
    public:
        // Texture Interface
        virtual T Evaluate(const SurfaceInteraction&) const = 0;
        virtual ~Texture() {}
    };

    Float Lanczos(Float, Float tau = 2);
    Float Noise(Float x, Float y = .5f, Float z = .5f);
    Float Noise(const Point3f& p);
    Float FBm(const Point3f& p, const Vector3f& dpdx, const Vector3f& dpdy,
        Float omega, int octaves);
    Float Turbulence(const Point3f& p, const Vector3f& dpdx, const Vector3f& dpdy,
        Float omega, int octaves);

}  // namespace pbrt



#endif