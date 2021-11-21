#ifndef MYH_H
#define MYH_H

#include "rtweekend.h"
#include "camera.h"
#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "moving_sphere.h"
#include "material.h"
#include <iostream>
#include "aarect.h"
#include "box.h"
#include "constant_medium.h"
#include "bvh.h"
#include "pdf.h"
#include <windows.h>
#include <thread>
#include "spectrum.h"

void vec3toarray3(vec3 a, double b[3])
{
    for (int i = 0; i < 3; i++)
    {
        b[i] = a.e[i];
    }
}

void array3tovec3(vec3 a, double b[3])
{
    for (int i = 0; i < 3; i++)
    {
        a.e[i] = b[i];
    }
}

hittable_list cornell_box() {
    hittable_list objects;

    auto red = make_shared<lambertian>(color(.65, .05, .05));
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    auto green = make_shared<lambertian>(color(.12, .45, .15));
    auto light = make_shared<diffuse_light>(color(15, 15, 15));

    hittable_list lights;
    lights.add(make_shared<xz_rect>(213, 343, 227, 332, 554, nullptr));
    lights.add(make_shared<sphere>(point3(190, 90, 190), 90, nullptr));
    //lights.add(make_shared<yz_rect>(0, 555, 0, 555, 0, nullptr));

    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(make_shared<flip_face>(make_shared<xz_rect>(213, 343, 227, 332, 554, light)));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

    /*shared_ptr<material> aluminum = make_shared<metal>(color(0.8, 0.85, 0.88), 0.0);
    shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), aluminum);
    box1 = make_shared<rotate_y>(box1, 15);
    box1 = make_shared<translate>(box1, vec3(265, 0, 295));
    objects.add(box1);*/
    shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), white);
    box1 = make_shared<rotate_y>(box1, 15);
    box1 = make_shared<translate>(box1, vec3(265, 0, 295));
    objects.add(box1);

    /*shared_ptr<hittable> box2 = make_shared<box>(point3(0, 0, 0), point3(165, 165, 165), white);
    box2 = make_shared<rotate_y>(box2, -18);
    box2 = make_shared<translate>(box2, vec3(130, 0, 65));
    objects.add(box2);*/

    auto glass = make_shared<dielectric>(1.5);
    objects.add(make_shared<sphere>(point3(190, 90, 190), 90, glass));

    return objects;
}

bool ifspectrum = false;
hittable_list cornell_box_spectrum() {

    ifspectrum = true;

    double white_array[3];
    double green_array[3];
    double red_array[3];
    double Light[3];
    vec3toarray3(color(.73, .73, .73), white_array);
    vec3toarray3(color(.65, .05, .05), red_array);
    vec3toarray3(color(.12, .45, .15), green_array);
    vec3toarray3(color(15, 15, 15), Light);
    /*std::cout << white_array[0] << "\n";
    std::cout << green_array[0] << "\n";
    std::cout << red_array[0] << "\n";
    std::cout << Light[0] << "\n";*/

    //pbrt::SampledSpectrum Lc = Lc.FromRGB(Light, pbrt::SpectrumType::Illuminant);
    pbrt::SampledSpectrum basic;
    basic.Init();
    pbrt::SampledSpectrum Lc = basic.FromRGB(Light, pbrt::SpectrumType::Illuminant);
    pbrt::SampledSpectrum Wc = basic.FromRGB(white_array, pbrt::SpectrumType::Reflectance);
    pbrt::SampledSpectrum Gc = basic.FromRGB(green_array, pbrt::SpectrumType::Reflectance);
    pbrt::SampledSpectrum Rc = basic.FromRGB(red_array, pbrt::SpectrumType::Reflectance);

    Lc.ToXYZ(Light);
    Wc.ToXYZ(white_array);
    Gc.ToXYZ(green_array);
    Rc.ToXYZ(red_array);
   /* std::cout << white_array[0] << "\n";
    std::cout << green_array[0] << "\n";
    std::cout << red_array[0] << "\n";
    std::cout << Light[0] << "\n";*/
    vec3 xyzwc(white_array[0], white_array[1], white_array[2]);
    vec3 xyzred(red_array[0], red_array[1], red_array[2]);
    vec3 xyzgre(green_array[0], green_array[1], green_array[2]);
    vec3 xyzL(Light[0], Light[1], Light[2]);

    /*array3tovec3(xyzwc, white_array);
    array3tovec3(xyzred, red_array);
    array3tovec3(xyzgre, green_array);
    array3tovec3(xyzL, Light);*/

    //std::cout << red_array[0];
    /*pbrt::RGBSpectrum whitecolor = whitecolor.FromRGB(white_array);
    pbrt::RGBSpectrum greencolor = greencolor.FromRGB(green_array);
    pbrt::RGBSpectrum redcolor = redcolor.FromRGB(red_array);*/
   
    /*std::cout << xyzwc << "\n";
    std::cout << xyzred << "\n";
    std::cout << xyzgre << "\n";
    std::cout << xyzL << "\n";*/



    hittable_list objects;

    auto red = make_shared<lambertian>(xyzred);
    auto white = make_shared<lambertian>(xyzwc);
    auto green = make_shared<lambertian>(xyzgre);
    auto light = make_shared<diffuse_light>(xyzL);

    hittable_list lights;
    lights.add(make_shared<xz_rect>(213, 343, 227, 332, 554, nullptr));
    lights.add(make_shared<sphere>(point3(190, 90, 190), 90, nullptr));

    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(make_shared<flip_face>(make_shared<xz_rect>(213, 343, 227, 332, 554, light)));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

    /*shared_ptr<material> aluminum = make_shared<metal>(color(0.8, 0.85, 0.88), 0.0);
    shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), aluminum);
    box1 = make_shared<rotate_y>(box1, 15);
    box1 = make_shared<translate>(box1, vec3(265, 0, 295));
    objects.add(box1);*/
    shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), white);
    box1 = make_shared<rotate_y>(box1, 15);
    box1 = make_shared<translate>(box1, vec3(265, 0, 295));
    objects.add(box1);

    /*shared_ptr<hittable> box2 = make_shared<box>(point3(0, 0, 0), point3(165, 165, 165), white);
    box2 = make_shared<rotate_y>(box2, -18);
    box2 = make_shared<translate>(box2, vec3(130, 0, 65));
    objects.add(box2);*/

    auto glass = make_shared<dielectric>(1.5);
    objects.add(make_shared<sphere>(point3(190, 90, 190), 90, glass));

    return objects;
}

//hittable_list cornell_smoke() {
//    hittable_list objects;
//
//    auto red = make_shared<lambertian>(color(.65, .05, .05));
//    auto white = make_shared<lambertian>(color(.73, .73, .73));
//    auto green = make_shared<lambertian>(color(.12, .45, .15));
//    auto light = make_shared<diffuse_light>(color(7, 7, 7));
//
//    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
//    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
//    objects.add(make_shared<xz_rect>(113, 443, 127, 432, 554, light));
//    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
//    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
//    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));
//
//    shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), white);
//    box1 = make_shared<rotate_y>(box1, 15);
//    box1 = make_shared<translate>(box1, vec3(265, 0, 295));
//
//    shared_ptr<hittable> box2 = make_shared<box>(point3(0, 0, 0), point3(165, 165, 165), white);
//    box2 = make_shared<rotate_y>(box2, -18);
//    box2 = make_shared<translate>(box2, vec3(130, 0, 65));
//
//    objects.add(make_shared<constant_medium>(box1, 0.01, color(0, 0, 0)));
//    objects.add(make_shared<constant_medium>(box2, 0.01, color(1, 1, 1)));
//
//    return objects;
//}

//hittable_list final_scene() {
//    hittable_list boxes1;
//    auto ground = make_shared<lambertian>(color(0.48, 0.83, 0.53));
//
//    const int boxes_per_side = 20;
//    for (int i = 0; i < boxes_per_side; i++) {
//        for (int j = 0; j < boxes_per_side; j++) {
//            auto w = 100.0;
//            auto x0 = -1000.0 + i * w;
//            auto z0 = -1000.0 + j * w;
//            auto y0 = 0.0;
//            auto x1 = x0 + w;
//            auto y1 = random_double(1, 101);
//            auto z1 = z0 + w;
//
//            boxes1.add(make_shared<box>(point3(x0, y0, z0), point3(x1, y1, z1), ground));
//        }
//    }
//
//    hittable_list objects;
//
//    //objects.add(make_shared<bvh_node>(boxes1, 0, 1));
//
//    auto light = make_shared<diffuse_light>(color(7, 7, 7));
//    objects.add(make_shared<xz_rect>(123, 423, 147, 412, 554, light));
//
//    hittable_list lights;
//    lights.add(make_shared<xz_rect>(123, 423, 147, 412, 554, nullptr));//
//
//
//    auto center1 = point3(400, 400, 200);
//    auto center2 = center1 + vec3(30, 0, 0);
//    auto moving_sphere_material = make_shared<lambertian>(color(0.7, 0.3, 0.1));
//    //objects.add(make_shared<moving_sphere>(center1, center2, 0, 1, 50, moving_sphere_material));
//
//    objects.add(make_shared<sphere>(point3(260, 150, 45), 50, make_shared<dielectric>(1.5)));
//    objects.add(make_shared<sphere>(
//        point3(0, 150, 145), 50, make_shared<metal>(color(0.8, 0.8, 0.9), 1.0)
//        ));
//
//    auto boundary = make_shared<sphere>(point3(360, 150, 145), 70, make_shared<dielectric>(1.5));
//    objects.add(boundary);
//    //objects.add(make_shared<constant_medium>(boundary, 0.2, color(0.2, 0.4, 0.9)));
//    boundary = make_shared<sphere>(point3(0, 0, 0), 5000, make_shared<dielectric>(1.5));
//    //objects.add(make_shared<constant_medium>(boundary, .0001, color(1, 1, 1)));
//
//    auto emat = make_shared<lambertian>(make_shared<image_texture>("earthmap.jpg"));
//    objects.add(make_shared<sphere>(point3(400, 200, 400), 100, emat));
//    auto pertext = make_shared<noise_texture>(0.1);
//    //objects.add(make_shared<sphere>(point3(220, 280, 300), 80, make_shared<lambertian>(pertext)));
//
//    hittable_list boxes2;
//    auto white = make_shared<lambertian>(color(.73, .73, .73));
//    int ns = 1000;
//    for (int j = 0; j < ns; j++) {
//        boxes2.add(make_shared<sphere>(point3::random(0, 165), 10, white));
//    }
//
//    objects.add(make_shared<translate>(
//        make_shared<rotate_y>(
//            make_shared<bvh_node>(boxes2, 0.0, 1.0), 15),
//        vec3(-100, 270, 395)
//        )
//    );
//
//    return objects;
//}



color ray_color(const ray& r, const color& background, const hittable& world,
    shared_ptr<hittable>& lights, int depth) {
    hit_record rec;
    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);

    //这里0.001能极大提速并且得到一致的效果，原因是在碰撞点射出的光线在t=0时刻总是和物体相交？
// If the ray hits nothing, return the background color.
    if (!world.hit(r, 0.001, infinity, rec))
        return background;

    scatter_record srec;
    color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
    if (!rec.mat_ptr->scatter(r, rec, srec))
        return emitted;

    if (srec.is_specular) {
        return srec.attenuation
            * ray_color(srec.specular_ray, background, world, lights, depth - 1);
    }

    auto light_ptr = make_shared<hittable_pdf>(lights, rec.p);
    mixture_pdf p(light_ptr, srec.pdf_ptr);

    ray scattered = ray(rec.p, p.generate(), r.time());
    auto pdf_val = p.value(scattered.direction());

    return emitted
        + srec.attenuation * rec.mat_ptr->scattering_pdf(r, rec, scattered)
        * ray_color(scattered, background, world, lights, depth - 1) / pdf_val;
}


// Image


const auto aspect_ratio = 1.0 / 1.0;
const int image_width = 600;
const int image_height = static_cast<int>(image_width / aspect_ratio);
const int samples_per_pixel = 100;
const int max_depth = 50;

//world

//康奈尔盒子
auto world = cornell_box_spectrum();
//auto world = cornell_box();
color background(0, 0, 0);
shared_ptr<hittable> lights =
make_shared<xz_rect>(213, 343, 227, 332, 554, shared_ptr<material>());
//康奈尔盒子

//hittable_list lights;
//lights.add(make_shared<xz_rect>(213, 343, 227, 332, 554, 0));
//lights.add(make_shared<sphere>(point3(190, 90, 190), 90, 0));
//auto world = final_scene();
//shared_ptr<hittable> lights =
//make_shared<xz_rect>(123, 423, 147, 412, 554, shared_ptr<material>());
//color background(0, 0, 0);

/*  switch (0) {
  case 1:
      world = random_scene();
      background = color(0.70, 0.80, 1.00);
      lookfrom = point3(13, 2, 3);
      lookat = point3(0, 0, 0);
      vfov = 20.0;
      aperture = 0.1;
      break;

  case 2:
      world = two_spheres();
      background = color(0.70, 0.80, 1.00);
      lookfrom = point3(13, 2, 3);
      lookat = point3(0, 0, 0);
      vfov = 20.0;
      break;

  case 3:
      world = two_perlin_spheres();
      background = color(0.70, 0.80, 1.00);
      lookfrom = point3(13, 2, 3);
      lookat = point3(0, 0, 0);
      vfov = 20.0;
      break;

  case 4:
      world = earth();
      background = color(0.70, 0.80, 1.00);
      lookfrom = point3(13, 2, 3);
      lookat = point3(0, 0, 0);
      vfov = 20.0;
      break;

  case 5:
      world = sphere_light();
      samples_per_pixel = 400;
      background = color(0, 0, 0);
      lookfrom = point3(26, 3, 6);
      lookat = point3(0, 2, 0);
      vfov = 20.0;
      break;

  case 6:
      world = cornell_box();
      aspect_ratio = 1.0;
      image_width = 600;
      samples_per_pixel = 200;
      background = color(0, 0, 0);
      lookfrom = point3(278, 278, -800);
      lookat = point3(278, 278, 0);
      vfov = 40.0;
      break;

  case 7:
      world = cornell_smoke();
      aspect_ratio = 1.0;
      image_width = 600;
      samples_per_pixel = 200;
      lookfrom = point3(278, 278, -800);
      lookat = point3(278, 278, 0);
      vfov = 40.0;
      break;

  default:
  case 8:
      world = final_scene();
      aspect_ratio = 1.0;
      image_width = 800;
      samples_per_pixel = 10000;
      background = color(0, 0, 0);
      lookfrom = point3(478, 278, -600);
      lookat = point3(278, 278, 0);
      vfov = 40.0;
      break;
  }*/

  // Camera
point3 lookfrom(278, 278, -800);
point3 lookat(278, 278, 0);


vec3 vup(0, 1, 0);
auto dist_to_focus = 10.0;
auto aperture = 0.0;
auto vfov = 40.0;
auto time0 = 0.0;
auto time1 = 1.0;

camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, time0, time1); 

//
std::vector<std::vector<vec3>> picture;
int dt = image_height / 4;
int t0 = dt;
int t1 = 2 * dt;
int t2 = 3 * dt;


#endif
