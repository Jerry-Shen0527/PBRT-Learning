#ifndef TEXTURE_H
#define TEXTURE_H

#include "rtweekend.h"
#include "rtw_stb_image.h"
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

class texture {
public:
    virtual color value(double u, double v, const point3& p) const = 0;
};

class solid_color : public texture {
public:
    solid_color() {}
    solid_color(color c) : color_value(c) {}

    solid_color(double red, double green, double blue)
        : solid_color(color(red, green, blue)) {}

    virtual color value(double u, double v, const vec3& p) const override {
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
            : even(make_shared<solid_color>(c1)) , odd(make_shared<solid_color>(c2)) {}

        virtual color value(double u, double v, const point3& p) const override {
            auto sines = sin(10*p.x())*sin(10*p.y())*sin(10*p.z());
            //auto sines = sin(100 * u) * sin(100 * v);
            if (sines < 0)
                return odd->value(u, v, p);
            else
                return even->value(u, v, p);
        }

    public:
        shared_ptr<texture> odd;
        shared_ptr<texture> even;
};

class perlin {
public:
    perlin() {
        ranvec = new vec3[point_count];
        for (int i = 0; i < point_count; ++i) {
            ranvec[i] = unit_vector(vec3::random(-1, 1));
        }


        perm_x = perlin_generate_perm();
        perm_y = perlin_generate_perm();
        perm_z = perlin_generate_perm();
    }

    ~perlin() {
        delete[] ranvec;
        delete[] perm_x;
        delete[] perm_y;
        delete[] perm_z;
    }

    double noise(const point3& p) const {
        auto u = p.x() - floor(p.x());
        auto v = p.y() - floor(p.y());
        auto w = p.z() - floor(p.z());
        auto i = static_cast<int>(floor(p.x()));
        auto j = static_cast<int>(floor(p.y()));
        auto k = static_cast<int>(floor(p.z()));
        vec3 c[2][2][2];

        for (int di = 0; di < 2; di++)
            for (int dj = 0; dj < 2; dj++)
                for (int dk = 0; dk < 2; dk++)
                    c[di][dj][dk] = ranvec[
                        perm_x[(i + di) & 255] ^
                            perm_y[(j + dj) & 255] ^
                            perm_z[(k + dk) & 255]
                    ];

        return perlin_interp(c, u, v, w);
    }

    double turb(const point3& p, int depth = 7) const {
        auto accum = 0.0;
        auto temp_p = p;
        auto weight = 1.0;

        for (int i = 0; i < depth; i++) {
            accum += weight * noise(temp_p);
            weight *= 0.5;
            temp_p *= 2;
        }

        return fabs(accum);
    }

private:
    static const int point_count = 256;
    double* ranfloat;
    int* perm_x;
    int* perm_y;
    int* perm_z;
    vec3* ranvec;

    static int* perlin_generate_perm() {
        auto p = new int[point_count];

        for (int i = 0; i < perlin::point_count; i++)
            p[i] = i;

        permute(p, point_count);

        return p;
    }

    static void permute(int* p, int n) {
        for (int i = n - 1; i > 0; i--) {
            int target = random_int(0, i);
            int tmp = p[i];
            p[i] = p[target];
            p[target] = tmp;
        }
    }

    static double trilinear_interp(double c[2][2][2], double u, double v, double w) {
        auto accum = 0.0;
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++)
                    accum += (i * u + (1 - i) * (1 - u)) *
                    (j * v + (1 - j) * (1 - v)) *
                    (k * w + (1 - k) * (1 - w)) * c[i][j][k];

        return accum;
    }

    static double perlin_interp(vec3 c[2][2][2], double u, double v, double w) {
        auto uu = u * u * (3 - 2 * u);
        auto vv = v * v * (3 - 2 * v);
        auto ww = w * w * (3 - 2 * w);
        auto accum = 0.0;

        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++) {
                    vec3 weight_v(u - i, v - j, w - k);
                    accum += (i * uu + (1 - i) * (1 - uu))
                        * (j * vv + (1 - j) * (1 - vv))
                        * (k * ww + (1 - k) * (1 - ww))
                        * dot(c[i][j][k], weight_v);
                }

        return accum;
    }
};

class noise_texture : public texture {
public:
    noise_texture() {}
    noise_texture(double sc) : scale(sc) {}

    virtual color value(double u, double v, const point3& p) const override {
        return color(1, 1, 1) * 0.5 * (1 + sin(scale * p.z() + 10 * noise.turb(p)));
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

    image_texture(cv::String filename) {
     
        image_mat = cv::imread(filename);
        cv::cvtColor(image_mat, image_mat, cv::COLOR_BGR2RGB);
        width = image_mat.cols;
        height = image_mat.rows;
       /* auto components_per_pixel = bytes_per_pixel;

        data = stbi_load(
            filename, &width, &height, &components_per_pixel, components_per_pixel);

        if (!data) {
            std::cerr << "ERROR: Could not load texture image file '" << filename << "'.\n";
            width = height = 0;
        }

        bytes_per_scanline = bytes_per_pixel * width;*/
    }

    ~image_texture() {
        delete data;
    }
    virtual color value(double u, double v, const vec3& p) const override {
        // Clamp input texture coordinates to [0,1] x [1,0]
        u = clamp(u, 0.0, 1.0);
        v = 1.0 - clamp(v, 0.0, 1.0);  // Flip V to image coordinates

        auto i = static_cast<int>(u * width);
        auto j = static_cast<int>(v * height);

        // Clamp integer mapping, since actual coordinates should be less than 1.0
        if (i >= width)  i = width - 1;
        if (j >= height) j = height - 1;
        const auto color_scale = 1.0 / 255.0;
       // std::cout << image_mat.at<cv::Vec3b>(j, i) << std::endl;
        return color(color_scale * image_mat.at<cv::Vec3b>(j, i)[0], 
                     color_scale * image_mat.at<cv::Vec3b>(j, i)[1], 
                     color_scale * image_mat.at<cv::Vec3b>(j, i)[2]);
    }
/*
    virtual color value(double u, double v, const vec3& p) const override {
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

        return color(color_scale * pixel[0], color_scale * pixel[1], color_scale * pixel[2]);
    }*/

private:
    cv::Mat image_mat;
    unsigned char* data;
    int width, height;
    int bytes_per_scanline;
};

#endif