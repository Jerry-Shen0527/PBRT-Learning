#include "vec3.h"
#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "bvh.h"

#include<iostream>
#include<thread>
#include<vector>
#include<ctime>
#include "omp.h"
using namespace std;

//#include "rtweekend.h"
#include "camera.h"
#include "material.h"
#include "moving_sphere.h"
#include "aarect.h"
#include "box.h"
#include "constant_medium.h"
#include "pdf.h"
//#include "my_image.h"
#include "spectrum.h"
//#include "geometry.h"
#include"transform.h"

color ray_color(
    const ray& r, const color& background, const hittable& world,
    const shared_ptr<hittable>& lights, int depth
) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);

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

hittable_list random_scene() {
    hittable_list world;

    auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));
    world.add(make_shared<sphere>(Point3f(0, -1000, 0), 1000, make_shared<lambertian>(checker)));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_Float();
            Point3f center(a + 0.9 * random_Float(), 0.2, b + 0.9 * random_Float());

            if ((center - Point3f(4, 0.2, 0)).Length() > 0.9) {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    auto center2 = center + Vector3f(0, random_Float(0, .5), 0);
                    world.add(make_shared<moving_sphere>(
                        center, center2, 0.0, 1.0, 0.2, sphere_material));
                }
                else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = random_Float(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(Point3f(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(Point3f(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(Point3f(4, 1, 0), 1.0, material3));

    return world;
}

hittable_list two_spheres() {
    hittable_list objects;

    auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));

    objects.add(make_shared<sphere>(Point3f(0, -10, 0), 10, make_shared<lambertian>(checker)));
    objects.add(make_shared<sphere>(Point3f(0, 10, 0), 10, make_shared<lambertian>(checker)));

    return objects;
}

hittable_list two_perlin_spheres() {
    hittable_list objects;

    auto pertext = make_shared<noise_texture>(4);
    objects.add(make_shared<sphere>(Point3f(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
    objects.add(make_shared<sphere>(Point3f(0, 2, 0), 2, make_shared<lambertian>(pertext)));

    return objects;
}

hittable_list earth() {
    auto earth_texture = make_shared<image_texture>("earthmap.jpg");
    auto earth_surface = make_shared<lambertian>(earth_texture);
    auto globe = make_shared<sphere>(Point3f(0, 0, 0), 2, earth_surface);

    return hittable_list(globe);
}

hittable_list simple_light() {
    hittable_list objects;

    auto pertext = make_shared<noise_texture>(4);
    objects.add(make_shared<sphere>(Point3f(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
    objects.add(make_shared<sphere>(Point3f(0, 2, 0), 2, make_shared<lambertian>(pertext)));

    auto difflight = make_shared<diffuse_light>(color(4, 4, 4));
    objects.add(make_shared<xy_rect>(3, 5, 1, 3, -2, difflight));
    objects.add(make_shared<sphere>(Point3f(0, 7, 2), 2, difflight));

    return objects;
}

hittable_list cornell_box() {
    hittable_list objects;

    auto red = make_shared<lambertian>(color(.65, .05, .05));
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    auto green = make_shared<lambertian>(color(.12, .45, .15));
    auto light = make_shared<diffuse_light>(color(15, 15, 15));

    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(make_shared<flip_face>(make_shared<xz_rect>(213, 343, 227, 332, 554, light)));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

    //shared_ptr<material> aluminum = make_shared<metal>(color(0.8, 0.85, 0.88), 0.0);
    //shared_ptr<hittable> box1 = make_shared<box>(Point3f(0, 0, 0), Point3f(165, 330, 165), aluminum);
    shared_ptr<hittable> box1 = make_shared<box>(Point3f(0, 0, 0), Point3f(165, 330, 165), white);
    box1 = make_shared<rotate_y>(box1, 15);
    box1 = make_shared<translate>(box1, Vector3f(265, 0, 295));
    objects.add(box1);

    auto glass = make_shared<dielectric>(1.5);
    objects.add(make_shared<sphere>(Point3f(190, 90, 190), 90, glass));

    return objects;
}

hittable_list cornell_box_spectrum() {
    Float xyzred[3] = {};
    Float xyzwhite[3] = {};
    Float xyzgreen[3] = {};
    Float xyzlight[3] = {};
    
    Float rgbred[3] = { .65, .05, .05 };
    Float rgbwhite[3] = { .73, .73, .73 };
    Float rgbgreen[3] = { .12, .45, .15 };
    Float rgblight[3] = { 15, 15, 15 };

    pbrt::SampledSpectrum::Init();
    pbrt::SampledSpectrum Samp;
    Samp = pbrt::SampledSpectrum::FromRGB(rgbred, pbrt::SpectrumType::Reflectance);
    Samp.ToXYZ(xyzred);
    Samp = pbrt::SampledSpectrum::FromRGB(rgbwhite, pbrt::SpectrumType::Reflectance);
    Samp.ToXYZ(xyzwhite);
    Samp = pbrt::SampledSpectrum::FromRGB(rgbgreen, pbrt::SpectrumType::Reflectance);
    Samp.ToXYZ(xyzgreen);
    Samp = pbrt::SampledSpectrum::FromRGB(rgblight, pbrt::SpectrumType::Illuminant);
    Samp.ToXYZ(xyzlight);

    hittable_list objects;

    auto xyz2vec = [](Float xyz[3]) {
        return Vector3f(xyz[0], xyz[1], xyz[2]);
    };

    auto red = make_shared<lambertian>(xyz2vec(xyzred));
    auto white = make_shared<lambertian>(xyz2vec(xyzwhite));
    auto green = make_shared<lambertian>(xyz2vec(xyzgreen));
    auto light = make_shared<diffuse_light>(xyz2vec(xyzlight));

    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(make_shared<flip_face>(make_shared<xz_rect>(213, 343, 227, 332, 554, light)));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

    //shared_ptr<material> aluminum = make_shared<metal>(color(0.8, 0.85, 0.88), 0.0);
    //shared_ptr<hittable> box1 = make_shared<box>(Point3f(0, 0, 0), Point3ff(165, 330, 165), aluminum);
    shared_ptr<hittable> box1 = make_shared<box>(Point3f(0, 0, 0), Point3f(165, 330, 165), white);
    box1 = make_shared<rotate_y>(box1, 15);
    box1 = make_shared<translate>(box1, Vector3f(265, 0, 295));
    objects.add(box1);

    auto glass = make_shared<dielectric>(1.5);
    objects.add(make_shared<sphere>(Point3f(190, 90, 190), 90, glass));

    return objects;
}

hittable_list cornell_smoke() {
    hittable_list objects;

    auto red = make_shared<lambertian>(color(.65, .05, .05));
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    auto green = make_shared<lambertian>(color(.12, .45, .15));
    auto light = make_shared<diffuse_light>(color(7, 7, 7));

    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(make_shared<xz_rect>(113, 443, 127, 432, 554, light));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

    shared_ptr<hittable> box1 = make_shared<box>(Point3f(0, 0, 0), Point3f(165, 330, 165), white);
    box1 = make_shared<rotate_y>(box1, 15);
    box1 = make_shared<translate>(box1, Vector3f(265, 0, 295));

    shared_ptr<hittable> box2 = make_shared<box>(Point3f(0, 0, 0), Point3f(165, 165, 165), white);
    box2 = make_shared<rotate_y>(box2, -18);
    box2 = make_shared<translate>(box2, Vector3f(130, 0, 65));

    objects.add(make_shared<constant_medium>(box1, 0.01, color(0, 0, 0)));
    objects.add(make_shared<constant_medium>(box2, 0.01, color(1, 1, 1)));

    return objects;
}

hittable_list final_scene() {
    hittable_list boxes1;
    auto ground = make_shared<lambertian>(color(0.48, 0.83, 0.53));

    const int boxes_per_side = 20;
    for (int i = 0; i < boxes_per_side; i++) {
        for (int j = 0; j < boxes_per_side; j++) {
            auto w = 100.0;
            auto x0 = -1000.0 + i * w;
            auto z0 = -1000.0 + j * w;
            auto y0 = 0.0;
            auto x1 = x0 + w;
            auto y1 = random_Float(1, 101);
            auto z1 = z0 + w;

            boxes1.add(make_shared<box>(Point3f(x0, y0, z0), Point3f(x1, y1, z1), ground));
        }
    }

    hittable_list objects;

    objects.add(make_shared<bvh_node>(boxes1, 0, 1));

    auto light = make_shared<diffuse_light>(color(7, 7, 7));
    objects.add(make_shared<xz_rect>(123, 423, 147, 412, 554, light));

    auto center1 = Point3f(400, 400, 200);
    auto center2 = center1 + Vector3f(30, 0, 0);
    auto moving_sphere_material = make_shared<lambertian>(color(0.7, 0.3, 0.1));
    objects.add(make_shared<moving_sphere>(center1, center2, 0, 1, 50, moving_sphere_material));

    objects.add(make_shared<sphere>(Point3f(260, 150, 45), 50, make_shared<dielectric>(1.5)));
    objects.add(make_shared<sphere>(
        Point3f(0, 150, 145), 50, make_shared<metal>(color(0.8, 0.8, 0.9), 1.0)
        ));

    auto boundary = make_shared<sphere>(Point3f(360, 150, 145), 70, make_shared<dielectric>(1.5));
    objects.add(boundary);
    objects.add(make_shared<constant_medium>(boundary, 0.2, color(0.2, 0.4, 0.9)));
    boundary = make_shared<sphere>(Point3f(0, 0, 0), 5000, make_shared<dielectric>(1.5));
    objects.add(make_shared<constant_medium>(boundary, .0001, color(1, 1, 1)));

    auto emat = make_shared<lambertian>(make_shared<image_texture>("earthmap.jpg"));
    objects.add(make_shared<sphere>(Point3f(400, 200, 400), 100, emat));
    auto pertext = make_shared<noise_texture>(0.1);
    objects.add(make_shared<sphere>(Point3f(220, 280, 300), 80, make_shared<lambertian>(pertext)));

    hittable_list boxes2;
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    int ns = 1000;
    for (int j = 0; j < ns; j++) {
        boxes2.add(make_shared<sphere>(Point3f::random(0, 165), 10, white));
    }

    objects.add(make_shared<translate>(
        make_shared<rotate_y>(
            make_shared<bvh_node>(boxes2, 0.0, 1.0), 15),
        Vector3f(-100, 270, 395)
        )
    );

    return objects;
}

int main() {

    // Image
    auto aspect_ratio = 16.0 / 9.0;
    int image_width = 400;
    //int image_height = static_cast<int>(image_width / aspect_ratio);
    int samples_per_pixel = 100;

    const int max_depth = 50;

    // World

    hittable_list world;

    Point3f lookfrom;
    Point3f lookat;
    auto vfov = 40.0;
    auto aperture = 0.0;
    color background(0, 0, 0);

    bool Isspectrum = false;
    switch (6) {
    case 1:
        world = random_scene();
        background = color(0.70, 0.80, 1.00);
        lookfrom = Point3f(13, 2, 3);
        lookat = Point3f(0, 0, 0);
        vfov = 20.0;
        aperture = 0.1;
        break;

    case 2:
        world = two_spheres();
        background = color(0.70, 0.80, 1.00);
        lookfrom = Point3f(13, 2, 3);
        lookat = Point3f(0, 0, 0);
        vfov = 20.0;
        break;

    case 3:
        world = two_perlin_spheres();
        background = color(0.70, 0.80, 1.00);
        lookfrom = Point3f(13, 2, 3);
        lookat = Point3f(0, 0, 0);
        vfov = 20.0;
        break;

    case 4:
        world = earth();
        background = color(0.70, 0.80, 1.00);
        lookfrom = Point3f(13, 2, 3);
        lookat = Point3f(0, 0, 0);
        vfov = 20.0;
        break;

    case 5:
        world = simple_light();
        samples_per_pixel = 400;
        background = color(0, 0, 0);
        lookfrom = Point3f(26, 3, 6);
        lookat = Point3f(0, 2, 0);
        vfov = 20.0;
        break;

    
    case 6:
        world = cornell_box();
        aspect_ratio = 1.0;
        image_width = 600;
        samples_per_pixel = 100;//200
        background = color(0, 0, 0);
        lookfrom = Point3f(278, 278, -800);
        lookat = Point3f(278, 278, 0);
        vfov = 40.0;
        break;

    case 7:
        world = cornell_smoke();
        aspect_ratio = 1.0;
        image_width = 600;
        samples_per_pixel = 200;
        lookfrom = Point3f(278, 278, -800);
        lookat = Point3f(278, 278, 0);
        vfov = 40.0;
        break;

    case 8:
        world = final_scene();
        aspect_ratio = 1.0;
        image_width = 800;
        //samples_per_pixel = 10000;
        samples_per_pixel = 800;
        background = color(0, 0, 0);
        lookfrom = Point3f(478, 278, -600);
        lookat = Point3f(278, 278, 0);
        vfov = 40.0;
        break;

    default:
    case 9:
        world = cornell_box_spectrum();
        aspect_ratio = 1.0;
        image_width = 600;
        samples_per_pixel = 100;//200
        background = color(0, 0, 0);
        lookfrom = Point3f(278, 278, -800);
        lookat = Point3f(278, 278, 0);
        vfov = 40.0;
        Isspectrum = true;
        break;
    }


    shared_ptr<hittable_list> lights = make_shared<hittable_list>();
    lights->add(make_shared<xz_rect>(213, 343, 227, 332, 554, shared_ptr<material>()));
    lights->add(make_shared<sphere>(Point3f(190, 90, 190), 90, shared_ptr<material>()));
    

    // Camera

    Vector3f vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    int image_height = static_cast<int>(image_width / aspect_ratio);

    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

    // Render
    //auto earth_texture = make_shared<my_image>("../PBRT-Learning/data/earthmap.jpg");
    //std::cout << earth_texture->get_width() << " " << earth_texture->get_height() << " " << earth_texture->get_bytes_per_scanline() << std::endl;
    
    const size_t core_num = std::thread::hardware_concurrency();
    size_t use_num = core_num;
    vector<thread> workers;
    workers.resize(use_num);
    
    thread ww;

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";
    
    size_t part_samples = samples_per_pixel / use_num;

    /*
    auto work = [=](int i, int j, size_t ps, color& pixel_color) {

        for (size_t s = 0; s < ps; ++s) {
            auto u = (i + random_Float()) / (image_width - 1);
            auto v = (j + random_Float()) / (image_height - 1);
            ray r = cam.get_ray(u, v);
            pixel_color += ray_color(r, background, world, lights, max_depth);
        }
    };*/
    
    //core code
    /*
    clock_t start, end;
    start = clock();
    for (int j = image_height - 1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
           
            //color pixel_color(0, 0, 0);
            //for (size_t idt = 0; idt < use_num; idt++)
             //   workers[idt] = thread(work, i, j, part_samples, std::ref(pixel_color));
            //for (auto& worker : workers)
             //   worker.join();

            color pixel_color(0, 0, 0);
#pragma omp parallel for
            for (int s = 0; s < samples_per_pixel; ++s) {
                   auto u = (i + random_Float()) / (image_width - 1);
                auto v = (j + random_Float()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, background, world, lights, max_depth);
            }
           
            if (Isspectrum)
            {
                Float xyz_pixel[3] = { pixel_color[0], pixel_color[1], pixel_color[2] };
                Float rgb_pixel[3] = {};
                pbrt::XYZToRGB(xyz_pixel, rgb_pixel);
                Vector3f rgb_color(rgb_pixel[0], rgb_pixel[1], rgb_pixel[2]);
                write_color(std::cout, rgb_color, samples_per_pixel);
            }
            else
                write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }
    end = clock();
    double run_time = (double)(end - start) / CLOCKS_PER_SEC;
    std::cerr << "\nrun_time: "<<run_time << " s.\nDone.\n";*/

    Point3f p1(1, 1, 1);
    cout <<"p1: " << p1 << endl;
    Vector3f vt(2.0, 3.0, 1.0);
    pbrt::Transform T = pbrt::Translate(vt);
    auto Rz = pbrt::RotateZ(90);
    auto Tp = T(p1);
    auto Rp = Rz(p1);
    Vector3f ve(0.1, 0.1, 0.1);
    auto pve = make_shared<Vector3f>(0.1, 0.1, 0.1);
    //auto vv = *pve;
    auto Tep = T(p1, pve.get());
    auto Rep = Rz(p1, &ve);
    cout << "TP: " << Tp << endl;
    cout << "TeP: " << Tep << endl;
    cout << "RP: " << Rp << endl;
    cout << "ReP: " << Rep << endl;
    auto RzT = Rz * T;
    auto RzTp = RzT(p1);
    cout << "PzTp: " << RzTp << endl;
    /*
    pbrt::SampledSpectrum::Init();
    color rgb1(0.5, 0.5, 0.5);
    pbrt::Float rgb[3] = { 0.75, 0.25, 0.8 };
    pbrt::Float torgb[3] = {};
    pbrt::Float xyz[3] = {};
    pbrt::Float xyztorgb[3] = {};
    pbrt::Float rgbtoxyz[3] = {};
    pbrt::SampledSpectrum sam = pbrt::SampledSpectrum::FromRGB(rgb, pbrt::SpectrumType::Reflectance);
   
    sam.ToXYZ(xyz);
    pbrt::RGBToXYZ(rgb, rgbtoxyz);

    sam.ToRGB(torgb);
    pbrt::XYZToRGB(rgbtoxyz, xyztorgb);

    cout << "original rgb \n";
    for (int i = 0; i < 3; i++) {
        cout << rgb[i] << " ";
    }
    cout << endl;

    cout << "sam from rgb to rgb \n";
    for (int i = 0; i < 3; i++) {
        cout << torgb[i] << " ";
    }
    cout << endl;

    cout << "rgb to xyz to rgb \n";
    for (int i = 0; i < 3; i++) {
        cout << xyztorgb[i] << " ";
    }
    cout << endl << "--------------" << endl;

    cout << "sam from rgb to xyz: " << endl;
    for (int i = 0; i < 3; i++) {
        cout << xyz[i] << " ";
    }

    cout << "\n rgbtoxyz: \n";
    for (int i = 0; i < 3; i++) {
        cout << rgbtoxyz[i] << " ";
    }
    cout << endl;
    */
}