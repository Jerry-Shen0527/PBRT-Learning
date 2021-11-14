#ifndef COLOR_H
#define COLOR_H

#include "vec3.h"

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>

void write_color(std::ostream &out, color pixel_color, int samples_per_pixel) {
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();
    // Replace NaN components with zero. See explanation in Ray Tracing: The Rest of Your Life.
    if (r != r) r = 0.0;
    if (g != g) g = 0.0;
    if (b != b) b = 0.0;
    // Divide the color by the number of samples and gamma-correct for gamma=2.0.
    auto scale = 1.0 / samples_per_pixel;
    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);

    // Write the translated [0,255] value of each color component.
    out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}

//use printf
void write_color(color pixel_color, int samples_per_pixel) {
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // Divide the color by the number of samples and gamma-correct for gamma=2.0.
    auto scale = 1.0 / samples_per_pixel;
    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);

    // Write the translated [0,255] value of each color component.
    printf("%d %d %d\n", 
        static_cast<int>(256 * clamp(r, 0.0, 0.999)), 
        static_cast<int>(256 * clamp(g, 0.0, 0.999)), 
        static_cast<int>(256 * clamp(b, 0.0, 0.999)));
}

void cv_write_color(cv::Mat& image_, int i, int j, const color& pixel_color, int samples_per_pixel)
{
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // Divide the color by the number of samples and gamma-correct for gamma=2.0.
    auto scale = 1.0 / samples_per_pixel;
    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);

    image_.at<cv::Vec3b>(j, i) = { static_cast<unsigned char>(256 * clamp(b, 0.0, 0.999)),
        static_cast<unsigned char>(256 * clamp(g, 0.0, 0.999)),
        static_cast<unsigned char>(256 * clamp(r, 0.0, 0.999)) };
}

inline double Gray(color c)
{
    return (c[0] + c[1] + c[2]) / 3;
}
#endif