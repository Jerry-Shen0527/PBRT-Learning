// PBRT.cpp: 定义应用程序的入口点。
//
#include <omp.h>
#include <time.h>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>

#include "rtweekend.h"
#include "color.h"
#include "hittable_list.h"
#include "sphere.h"
#include "triangle.h"
#include "camera.h"
#include "material.h"
#include "moving_sphere.h"
#include "texture.h"
#include "aarect.h"
#include "constant_medium.h"
#include "spectrum.h"
#include "transform.h"
#include "loadobj.h"
#include "primitive.h"

Color ray_color(const ray& r, const Color& background, const std::vector<shared_ptr<Primitive>> &obj, const hittable& world, shared_ptr<hittable>& lights) {
    if (r.depth >= 50)
        return Color(0.f);
    SurfaceInteraction isec;
    hit_record rec;

    bool is_isec = false;
    bool flag_obj = false;
    bool flag_world = false;
    for (int n = 0; n < obj.size(); n++)
    {
        if (obj[n]->Intersect(r, &isec))
        {
            flag_obj = true;
            is_isec = true;
        }
    }
    if (!world.hit(r, 0.001, r.tMax, rec))
    {   
        if (flag_obj == false)
            return background;
        else //intersect triangle obj
        {
            // light
            //  return Color(0.9f) /** ray_color(ray(r,true), background, obj, world, lights)*/;
            //lambert
            if(Dot(r.d,vec3(isec.n)) > 0)
                isec.n = -isec.n;
            
            auto light_ptr = make_shared<hittable_pdf>(lights, isec.p);
            auto scatter_pdf = make_shared<cosine_pdf>(vec3(isec.n));
            mixture_pdf p(light_ptr, scatter_pdf);

            ray scattered = ray(isec.p, p.generate(), r);
            auto pdf_val = p.value(scattered.direction());

            auto cosine = Dot(vec3(isec.n), unit_vector(scattered.direction()));
            cosine = cosine < 0 ? 0 : cosine / Pi;
            return Color::FromRGB(vec3(191,184,241)/255.0) * cosine
                * ray_color(scattered, background, obj, world, lights) / pdf_val;
        }
           
    }
    else
    {
        scatter_record srec;
        Color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
        if (!rec.mat_ptr->scatter(r, rec, srec))
            return emitted;
        if (srec.is_specular) {
            return srec.attenuation
                * ray_color(ray(srec.specular_ray, true), background, obj, world, lights);
        }
        auto light_ptr = make_shared<hittable_pdf>(lights, rec.p);
        mixture_pdf p(light_ptr, srec.pdf_ptr);

        ray scattered = ray(rec.p, p.generate(), r);
        auto pdf_val = p.value(scattered.direction());

        return emitted
            + srec.attenuation * rec.mat_ptr->scattering_pdf(r, rec, scattered)
            * ray_color(scattered, background, obj, world, lights) / pdf_val;

    }
    //if (flag_obj ==false)
    //    return Color(0.f);
    //else
    //    return Color::FromRGB(vec3(isec.n));
}

Color ray_color(const ray& r, const Color& background, const hittable& world,  shared_ptr<hittable>& lights) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (r.depth >= 50)
        return Color(0.f);

    // If the ray hits nothing, return the background color.
    if (!world.hit(r, 0.001, Infinity, rec))
        return background;
    scatter_record srec;
    Color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
    if (!rec.mat_ptr->scatter(r, rec, srec))
        return emitted;
    if (srec.is_specular) {
        return srec.attenuation
            * ray_color(ray(srec.specular_ray,true), background, world, lights);
    }
    auto light_ptr = make_shared<hittable_pdf>(lights, rec.p);
    mixture_pdf p(light_ptr, srec.pdf_ptr);

    ray scattered = ray(rec.p, p.generate(), r);
    auto pdf_val = p.value(scattered.direction());
   
    return emitted
        + srec.attenuation * rec.mat_ptr->scattering_pdf(r, rec, scattered)
        * ray_color(scattered, background, world, lights) / pdf_val;
}

extern std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(
    shared_ptr<Transform> ObjectToWorld, shared_ptr<Transform> WorldToObject,
    bool reverseOrientation, int nTriangles,
    const int* vertexIndices, int nVertices, const Point3f* p,
    const Vector3f* s, const Normal3f* n, const Point2f* uv,
    //const std::shared_ptr<Texture<Float>>& alphaMask),
    //const std::shared_ptr<Texture<Float>>& shadowAlphaMask,
    const int* faceIndices);

std::vector<shared_ptr<Primitive>> new_scene()
{
    shared_ptr<Transform> id = make_shared<Transform>();
    Sphere obj1(id, id, 1);
    auto S1 = GeometricPrimitive(make_shared<Sphere>(obj1));
    
    Transform Scale(Float x, Float y, Float z);
    
    shared_ptr<Transform> small = make_shared<Transform>(Scale(0.025, 0.025, 0.025)* Rotate(-90.0f, vec3(1,0,0)));
    //shared_ptr<Transform> big = make_shared<Transform>(Translate(vec3(365,100,200))*Scale(100, 100, 100) );
    Transform cube_pre = Scale(30, 30, 30);
    shared_ptr<Transform> big = make_shared<Transform>(Translate(vec3(-450,0,-100))*Scale(2,2,2)*Translate(vec3(365, 100, 200)) * Scale(2.5, 2.5, 2.5) * Rotate(-90.0f, vec3(1, 0, 0)) * Rotate(-90.0f, vec3(0, 0, 1)));
    shared_ptr<Transform> cube_trans = make_shared<Transform>((*big) * cube_pre* Rotate(45, vec3(0, 1, 0)));
    shared_ptr<Transform> hot_dog_trans = make_shared<Transform>(Translate(vec3(-450, 0, -100)) * Scale(2, 2, 2) * Translate(vec3(365, 100, 200)) * Scale(2.5, 2.5, 2.5) * Rotate(20,vec3(1,0,0))*Rotate(45, vec3(0, 0, 1)) * Rotate(180, vec3(0, 1, 0)));
    std::vector< shared_ptr<Primitive> > scene;
    Model qwq("D:\\QWQ\\data\\mesh\\triangle mesh\\Hot\ dog.obj");
    Mesh pwp = qwq.meshes[0];
    auto cube=CreateTriangleMesh(hot_dog_trans,make_shared<Transform>(Inverse(*hot_dog_trans)), false, pwp.f_num, pwp.f_indics, pwp.v_num, pwp.v_pos, pwp.vt, pwp.vn, pwp.uv, nullptr);
    for (auto& iter : cube)
    {
        scene.push_back(make_shared<GeometricPrimitive>(iter));
    }

    return { make_shared<pbrt_bvh_node>(scene, 0, scene.size() - 1) };
}

hittable_list test()
{
    hittable_list world;
    
    auto material3 = make_shared<lambertian>(color(0.8, 0.7, 0.7));
    world.add(make_shared<sphere>(point3(-4, -1, 0), 1.0, material3));
    //auto ground_t = make_shared<solid_color>(color(0.35, 0.3, 0.25));
   // world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(ground_t)));
    auto light = make_shared<diffuse_light>(color(4, 4, 4));
    world.add(make_shared<xz_rect>(-100,100, -100, 100, 10, light));

    return world;
}

hittable_list random_scene() {
    hittable_list world;

    auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));
    world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(checker)));

    int n = 11;
    for (int a = -n; a < n; a++) {
        for (int b = -n; b < n; b++) {
            auto choose_mat = RandomFloat();
            point3 center(a + 0.9 * RandomFloat(), 0.2, b + 0.9 * RandomFloat());

            if ((center - point3(4, 0.2, 0)).Length() > 0.9) {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    auto center2 = center + vec3(0, RandomFloat(0, .5), 0);
                    world.add(make_shared<moving_sphere>(
                        center, center2, 0.0, 1.0, 0.2, sphere_material));
                    world.setTime(0.0, 1.0);
                }
                else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = RandomFloat(0, 0.5);
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
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    return world;
}

hittable_list two_spheres() {
    hittable_list objects;
    color c1 = color(0.2, 0.3, 0.1);
    color c2 = color(17.0,64.0,148.0) / 255.0;
    color c3(0.8, 0, 0);
    auto checker = make_shared<checker_texture>(c3, color(0.9, 0.9, 0.9));
    //auto checker = make_shared<checker_texture>(color(0.3, 0.1, 0.1), color(0.9, 0.9, 0.9));
    
    objects.add(make_shared<sphere>(point3(0, -10, -5), 10, make_shared<lambertian>(checker)));
    objects.add(make_shared<sphere>(point3(0, 10, -5), 10, make_shared<lambertian>(checker)));

    return objects;
}

hittable_list two_perlin_spheres() {
    hittable_list objects;
    auto pertext = make_shared<noise_texture>(4);
    objects.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
    objects.add(make_shared<sphere>(point3(0, 2, 0), 2, make_shared<lambertian>(pertext)));
    objects.add(make_shared<sphere>(point3(0, 2, 0), 2, make_shared<lambertian>(pertext)));

    return objects;
}

hittable_list earth() {
    hittable_list objects;
    auto earth_texture = make_shared<image_texture>("E:\\PBRT\\PBRT-Learning\\image\\test_.png");
    auto earth_surface = make_shared<lambertian>(earth_texture);
    auto globe = make_shared<sphere>(point3(0, 0, 0), 2, earth_surface);
    objects.add(globe);
    auto difflight = make_shared<diffuse_light>(color(4, 4, 4));
    objects.add(make_shared<xy_rect>(3, 5, 1, 3, -2, difflight));
    return objects;
}

hittable_list cornell_box() {
    hittable_list objects;

    auto red = make_shared<lambertian>(color(.65, .05, .05));
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    auto green = make_shared<lambertian>(color(.12, .45, .15));
    auto light = make_shared<diffuse_light>(color(20, 18, 15)*2);

    //walls
    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(make_shared<xz_rect>(213, 343, 227, 332, 554, light));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

    ////blocks
    //shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), white);
    //Transform rotate1 = RotateX(15);
    //box1 = make_shared<rotate_y>(box1, 15);
    //box1 = make_shared<translate>(box1, vec3(265, 0, 295));
    //objects.add(box1);

    //shared_ptr<hittable> box2 = make_shared<box>(point3(0, 0, 0), point3(165, 165, 165), white);
    //box2 = make_shared<rotate_y>(box2, -18);
    //box2 = make_shared<translate>(box2, vec3(130, 0, 65));
    //objects.add(box2);

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

    shared_ptr<hittable> box1 = make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), white);
    box1 = make_shared<rotate_y>(box1, 15);
    box1 = make_shared<translate>(box1, vec3(265, 0, 295));

    shared_ptr<hittable> box2 = make_shared<box>(point3(0, 0, 0), point3(165, 165, 165), white);
    box2 = make_shared<rotate_y>(box2, -18);
    box2 = make_shared<translate>(box2, vec3(130, 0, 65));

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
            auto y1 = RandomFloat(1, 101);
            auto z1 = z0 + w;

            boxes1.add(make_shared<box>(point3(x0, y0, z0), point3(x1, y1, z1), ground));
        }
    }

    hittable_list objects;

    objects.add(make_shared<bvh_node>(boxes1, 0, 1));

    auto light = make_shared<diffuse_light>(color(7, 7, 7));
    objects.add(make_shared<xz_rect>(123, 423, 147, 412, 554, light));

    auto center1 = point3(400, 400, 200);
    auto center2 = center1 + vec3(30, 0, 0);
    auto moving_sphere_material = make_shared<lambertian>(color(0.7, 0.3, 0.1));
    objects.add(make_shared<moving_sphere>(center1, center2, 0, 1, 50, moving_sphere_material));

    objects.add(make_shared<sphere>(point3(260, 150, 45), 50, make_shared<dielectric>(1.5)));
    objects.add(make_shared<sphere>(
        point3(0, 150, 145), 50, make_shared<metal>(color(0.8, 0.8, 0.9), 1.0)
        ));

    auto boundary = make_shared<sphere>(point3(360, 150, 145), 70, make_shared<dielectric>(1.5));
    objects.add(boundary);
    objects.add(make_shared<constant_medium>(boundary, 0.2, color(0.2, 0.4, 0.9)));
    boundary = make_shared<sphere>(point3(0, 0, 0), 5000, make_shared<dielectric>(1.5));
    objects.add(make_shared<constant_medium>(boundary, .0001, color(1, 1, 1)));

    auto emat = make_shared<lambertian>(make_shared<image_texture>("E:\\PBRT\\PBRT-Learning\\image\\test_.png"));
    objects.add(make_shared<sphere>(point3(400, 200, 400), 100, emat));
   // auto pertext = make_shared<noise_texture>(0.1);
    //objects.add(make_shared<sphere>(point3(220, 280, 300), 80, make_shared<lambertian>(pertext)));

    hittable_list boxes2;
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    int ns = 1000;
    for (int j = 0; j < ns; j++) {
        boxes2.add(make_shared<sphere>(point3::random(0, 165), 10, white));
    }

    objects.add(make_shared<translate>(
        make_shared<rotate_y>(
            make_shared<bvh_node>(boxes2, 0.0, 1.0), 15),
        vec3(-100, 270, 395)
        )
    );

    return objects;
}


int main() {

    SampledSpectrum::Init();
    //std::cout << v << std::endl;
    // Image
    auto aspect_ratio = 16.0 / 9.0;
    int image_width = 400;
    int samples_per_pixel =20;
    const int max_depth = 50;

    // World

    hittable_list world;

    point3 lookfrom;
    point3 lookat;
    auto vfov = 40.0;
    auto aperture = 0.0;

    color background(0, 0, 0);

    switch (5) {
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
        samples_per_pixel = 100;
        background = color(0, 0, 0);
        lookfrom = point3(26, 3, 6);
        lookat = point3(0, 2, 0);
        vfov = 20.0;
        break;

    case 5:
        world = cornell_box();
        aspect_ratio = 1.0;
        image_width = 600;
        samples_per_pixel = 100;
        background = color(0, 0, 0);
        lookfrom = point3(278, 278, -800);
        lookat = point3(278, 278, 0);
        vfov = 40.0;
        break;

    case 6:
        world = cornell_smoke();
        aspect_ratio = 1.0;
        image_width = 600;
        samples_per_pixel = 200;
        lookfrom = point3(278, 278, -800);
        lookat = point3(278, 278, 0);
        vfov = 40.0;
        break;

    case 7:
        world = test();
        background = color(0.2, 0.2, 0.2);
        lookfrom = point3(13, 2, 3);
        lookat = point3(0, 0, 0);
        vfov = 20.0;
        break;

    default:
    case 8:
        world = final_scene();
        aspect_ratio = 1.0;
        image_width = 800;
        samples_per_pixel = 16;
        background = color(0, 0, 0);
        lookfrom = point3(478, 278, -600);
        lookat = point3(278, 278, 0);
        vfov = 40.0;
        break;
    }
    shared_ptr<hittable> lights = make_shared<xz_rect>(213, 343, 227, 332, 554, shared_ptr<material>());


    //bvh_node bvh(world, world._time0, world._time1);
    // Camera

    vec3 vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    int image_height = static_cast<int>(image_width / aspect_ratio);

    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

    // Render
    cv::Mat image_ = cv::Mat::zeros(image_height, image_width, CV_8UC3);
    clock_t start = clock();
    printf("P3\n%d %d\n255\n", image_width, image_height);
    Color background_sp = Color::FromRGB(background,SpectrumType::Illuminant);
    std::vector<shared_ptr<Primitive>> scene11 = new_scene();
    for (int j = image_height - 1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
#pragma omp parallel for
        for (int i = 0; i < image_width; ++i) {
            Color pixel_color(0.f);

            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + RandomFloat()) / (image_width - 1);
                auto v = (j + RandomFloat()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                //pixel_color += ray_color(r,background_sp, world,lights);
                pixel_color += ray_color(r, background_sp, scene11, world, lights);
            }
            //write_color( pixel_color, samples_per_pixel);
            cv_write_color(image_, i, image_height - 1 - j, pixel_color.ToColor(), samples_per_pixel);
        }
    }
    clock_t end = clock();
    std::cerr << "\nSpend time:" << (end - start) / 1000 << "s. Done.\n";

    cv::imwrite("E:\\PBRT\\PBRT-Learning\\image\\qwq.png", image_);
    //去噪
    cv::Mat result1, result2, result3, result4;
    blur(image_, result1, cv::Size(3, 3));
    cv::imwrite("E:\\PBRT\\PBRT-Learning\\image\\qwq_1.png", result1);

    GaussianBlur(image_, result2, cv::Size(3, 3), 0);
    cv::imwrite("E:\\PBRT\\PBRT-Learning\\image\\qwq_2.png", result2);

    medianBlur(image_, result3, 3);
    cv::imwrite("E:\\PBRT\\PBRT-Learning\\image\\qwq_3.png", result3);

    fastNlMeansDenoisingColored(image_, result4);
    cv::imwrite("E:\\PBRT\\PBRT-Learning\\image\\qwq_4.png", result4);
    
}


