//#include "rtweekend.h"
//#include "camera.h"
//#include "color.h"
//#include "hittable_list.h"
//#include "sphere.h"
//#include "moving_sphere.h"
//#include "material.h"
//#include <iostream>
//#include "aarect.h"
//#include "box.h"
//#include "constant_medium.h"
//#include "bvh.h"
//#include "pdf.h"
//#include <windows.h>
#include "myh.h"





void fun0()
{
    for (int j = t0 - 1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, background, world, lights, max_depth);
            }
            picture[j][i] = pixel_color;
        }
    }
}

void fun1()
{
    for (int j = t1 - 1; j >= t0; --j) {
        std::cerr << "\rScanlines remaining: " << j-t0 << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, background, world, lights, max_depth);
            }
            picture[j][i] = pixel_color;
           
        }
    }
}
void fun2()
{
    for (int j = t2 - 1; j >= t1; --j) {
        std::cerr << "\rScanlines remaining: " << j-t1 << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, background, world, lights, max_depth);
            }
            picture[j][i] = pixel_color;

        }
    }
}

void fun3()
{
    for (int j = image_height-1; j >= t2; --j) {
        std::cerr << "\rScanlines remaining: " << j -t2<< ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, background, world, lights, max_depth);
            }
            picture[j][i] = pixel_color;
            //std::cout << picture[j][i];
        }
    }
}

int main() {


    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    
    picture.resize(4*image_height);
    for (int j = 0; j < image_height; j++)
    {
        picture[j].resize(image_width * 2);
    }

    
    // Render
    //¶àÏß³Ì
    std::thread t1(fun0);
    std::thread t2(fun1);
    std::thread t3(fun2);
    std::thread t4(fun3);

    t1.join();
    t2.join();
    t3.join();
    t4.join();

    
    if (ifspectrum) {
        for (int j = image_height - 1; j >= 0; --j) {
            for (int i = 0; i < image_width; ++i) {
                double rgb[3];
                double xyz[3];
                vec3toarray3(picture[j][i], xyz);
                pbrt::XYZToRGB(xyz, rgb);
                write_color(std::cout, rgb, samples_per_pixel);
            }
        }
    }

    if (!ifspectrum) {
        for (int j = image_height - 1; j >= 0; --j) {
            for (int i = 0; i < image_width; ++i) {
                write_color(std::cout, picture[j][i], samples_per_pixel);
            }
        }
    }

    std::cerr << "\nDone.\n";
    
}
