#include "RaytracingPart3.h"

RaytracingPart3 ::RaytracingPart3(QWidget *parent)
    : QWidget(parent)
{
    ui.setupUi(this);
}

void RaytracingPart3::getImage() {
    filename = "../image/test5.jpg";

    //Image

    auto aspect_ratio = 16.0 / 9.0;
    int image_width = 400;
    int samples_per_pixel;
    const int max_depth = 2;

    //world
    Primitive_list world;

    Point3f lookfrom;
    Point3f lookat;
    auto vfov = 40.0;
    auto aperture = 0.0;
    color background(0, 0, 0);
    //Spectrum background(0.0);
    /*shared_ptr<Primitive> lights =
        //  make_shared<xz_rect>(213, 343, 227, 332, 554, shared_ptr<material>());
        make_shared<sphere>(Point3f(190, 90, 190), 90, shared_ptr<material>());*/

    auto lights = make_shared<Primitive_list>();
    auto red = make_shared<lambertian>(color(.65, .05, .05));
    lights->add(make_shared<xz_rect>(213, 343, 227, 332, 554, shared_ptr<material>()));

    /*std::shared_ptr<Primitive> accelerator =
        CreateKdTreeAccelerator(std::move((shared_ptr<Primitive>)lights));*/
    //lights->add(make_shared<sphere>(Point3f(190, 90, 190), 90, shared_ptr<material>()));

    /*shared_ptr<Primitive> lights =
        make_shared<xz_rect>(213, 343, 227, 332, 554, shared_ptr<material>());*/

    world = Test();
    

    aspect_ratio = 1.0;
    image_width = 600;
    samples_per_pixel = 100;
    background = color(1.0, 1.0, 1.0);
    lookfrom = Point3f(278, 278, -800);
    lookat = Point3f(278, 278, 0);
    vfov = 40.0;

    // Camera

    Vector3f vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    int image_height = static_cast<int>(image_width / aspect_ratio);

    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

    QImage image(image_width, image_height, QImage::Format_ARGB32);

    clock_t begin, duration;
    begin = clock();

#pragma omp parallel for
    for (int j = image_height - 1; j >= 0; --j) {
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);

            for (int s = 0; s < samples_per_pixel; ++s) {

                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                //cout << r.direction().x << " " << r.direction().y << " " << r.direction().z << endl;
                color c_tmp=ray_color(r, background, world, (shared_ptr<Primitive>)lights, max_depth);
               //color c_tmp = ray_color(r, background, world, max_depth);
                pixel_color += c_tmp;
                //cout << c_tmp.x << " " << c_tmp.y << " " << c_tmp.z << endl;

            }

            vector<int> rgb = write_color(pixel_color, samples_per_pixel);
            image.setPixel(QPoint(i, image_height - 1 - j), qRgb(rgb[0], rgb[1], rgb[2]));
        }
    }
    duration = clock() - begin;

    image.save(filename);
    cout << "Done:"<< (double)duration/ CLOCKS_PER_SEC << endl;
}


void RaytracingPart3::paintEvent(QPaintEvent*) {
    if (cal_bool) {
        //getImage();
        //getImageChecker();
        getImage();
        cal_bool = false;
    }

    QPainter p(this);
    p.translate(0, 0);
    p.drawImage(0, 0, QImage(filename));
}

double RaytracingPart3::hit_sphere(const Point3f& center, double radius, const ray& r) {
    Vector3f oc = r.origin() - center;
    auto a = r.direction().LengthSquared();
    auto half_b = Dot(oc, r.direction());
    auto c = oc.LengthSquared() - radius * radius;
    auto discriminant = half_b * half_b - a * c;

    if (discriminant < 0) {
        return -1.0;
    }
    else {
        return (-half_b - sqrt(discriminant)) / a;
    }
}


color RaytracingPart3::ray_color(
    const ray& r, const color& background, const Primitive& world,
    shared_ptr<Primitive>& lights, int depth) {
   
    SurfaceInteraction rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);

    // If the ray hits nothing, return the background color.
    if (!world.Intersect(r, 0.001, infinity, rec))
        return background;
    //cout << rec.p.x << " " << rec.p.y << " " << rec.p.z << endl;
    scatter_record srec;
    color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
    if (!rec.mat_ptr->scatter(r, rec, srec))
        return emitted;

    if (srec.is_specular) {
        return srec.attenuation.MutiXYZ(ray_color(srec.specular_ray, background, world, lights, depth - 1));
            
    }
    //cout << rec.p.x << " " << rec.p.y << " " << rec.p.z << endl;
    auto light_ptr = make_shared<Primitive_pdf>(lights, rec.p);
    mixture_pdf p(light_ptr, srec.pdf_ptr);
    
   // mixture_pdf p(light_ptr, srec.pdf_ptr);
    ray scattered = ray(rec.p, p.generate(), r.time());
    auto pdf_val = p.value(scattered.direction());
    
    return emitted
        + rec.mat_ptr->scattering_pdf(r, rec, scattered)* srec.attenuation.MutiXYZ(ray_color(scattered, background, world, lights, depth - 1))
        / pdf_val;
}

color RaytracingPart3::ray_color(
    const ray& r, const Spectrum& background, const Primitive& world,
    shared_ptr<Primitive>& lights, int depth) {
    SurfaceInteraction rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);

    // If the ray hits nothing, return the background color.
    if (!world.Intersect(r, 0.001, infinity, rec))
        return background.ToColor();

    scatter_record srec;
    color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
    if (!rec.mat_ptr->scatter(r, rec, srec))
        return emitted;

    if (srec.is_specular) {
        return srec.attenuationSpe.ToColor().MutiXYZ(ray_color(srec.specular_ray, background, world, lights, depth - 1))
            ;
    }

    auto light_ptr = make_shared<Primitive_pdf>(lights, rec.p);
    mixture_pdf p(light_ptr, srec.pdf_ptr);

    ray scattered = ray(rec.p, p.generate(), r.time());
    auto pdf_val = p.value(scattered.direction());

    return emitted
        + rec.mat_ptr->scattering_pdf(r, rec, scattered)*(srec.attenuationSpe.ToColor()).MutiXYZ(ray_color(scattered, background, world, lights, depth - 1))
         / pdf_val;
}

color RaytracingPart3::ray_color(
    const ray& r, const color& background, const Primitive& world,
    int depth) {
    SurfaceInteraction rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);

    // If the ray hits nothing, return the background color.
    if (!world.Intersect(r, 0.001, infinity, rec))
        return background;
    //cout << rec.p.x << " " << rec.p.y << " " << rec.p.z << endl;
    scatter_record srec;
    color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
    if (!rec.mat_ptr->scatter(r, rec, srec))
        return emitted;

    if (srec.is_specular) {
        return srec.attenuation.MutiXYZ(ray_color(srec.specular_ray, background, world, depth - 1));

    }
    //cout << rec.p.x << " " << rec.p.y << " " << rec.p.z << endl;
    mixture_pdf p(srec.pdf_ptr, srec.pdf_ptr);
    
    ray scattered = ray(rec.p, p.generate(), r.time());
    auto pdf_val = p.value(scattered.direction());
    //cout << scattered.direction()[0] << " " << scattered.direction()[1] << " " << scattered.direction()[2] << endl;
    //cout << pdf_val <<" "<< light_ptr->value(scattered.direction()) << " " << srec.pdf_ptr->value(scattered.direction()) << endl;
    return emitted
        + rec.mat_ptr->scattering_pdf(r, rec, scattered) * srec.attenuation.MutiXYZ(ray_color(scattered, background, world, depth - 1))
        / pdf_val;
}

vector<int> RaytracingPart3::write_color(color pixel_color, int samples_per_pixel) {
    auto r = pixel_color.x;
    auto g = pixel_color.y;
    auto b = pixel_color.z;

    // Replace NaN components with zero. See explanation in Ray Tracing: The Rest of Your Life.
    if (r != r) r = 0.0;
    if (g != g) g = 0.0;
    if (b != b) b = 0.0;

    auto scale = 1.0 / samples_per_pixel;
    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);


    vector<int> rgb;
    // Write the translated [0,255] value of each color component.
    rgb.push_back(static_cast<int>(256 * clamp(r, 0.0, 0.999)));
    rgb.push_back(static_cast<int>(256 * clamp(g, 0.0, 0.999)));
    rgb.push_back(static_cast<int>(256 * clamp(b, 0.0, 0.999)));

    return rgb;
}


/*Primitive_list RaytracingPart3::cornell_box() {
    Primitive_list objects;

    auto red = make_shared<lambertian>(color(.65, .05, .05));
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    auto green = make_shared<lambertian>(color(.12, .45, .15));
    auto light = make_shared<diffuse_light>(color(10, 10, 10));

    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(make_shared<flip_face>(make_shared<xz_rect>(213, 343, 227, 332, 554, light)));
    //objects.add(make_shared<xz_rect>(213, 343, 227, 332, 554, light));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

    auto glass = make_shared<dielectric>(1.5);
    objects.add(make_shared<sphere>(Point3f(190, 90, 190), 80, green));

    std::shared_ptr<Primitive> accelerator =
        CreateBVHAccelerator(std::move(objects.objects));
    return accelerator;
}*/

Primitive_list RaytracingPart3::Test() {
    Primitive_list objects;

    auto red = make_shared<lambertian>(color(.65, .05, .05));
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    auto green = make_shared<lambertian>(color(.12, .45, .15));
    auto light = make_shared<diffuse_light>(color(10, 10, 10));
    auto yellow= make_shared<lambertian>(color(.12, .65, .65));

    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<flip_face>(make_shared<xz_rect>(213, 343, 227, 332, 554, light)));
   /* objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<flip_face>(make_shared<xz_rect>(213, 343, 227, 332, 554, red)));
    objects.add(make_shared<xz_rect>(213, 343, 227, 332, 554, red));*/

    //objects.add(make_shared<sphere>(Point3f(190, 90, 190), 90, red));
    /*Transform t = Translate(Point3f(190, 90, 190));
    //cout << t.GetMatrix().m[0][3] << endl;
    Transform rt = Translate(Point3f(-190, -90, -190));

    vector<Point3f> pos;
    vector<int> face;
    vector<shared_ptr<material>> mat_ptr;
    mat_ptr.push_back(red);
    mat_ptr.push_back(green);
    pos.push_back(Point3f(100, 70, 160));
    pos.push_back(Point3f(210, 150, 210));
    pos.push_back(Point3f(210, 70, 70));
    pos.push_back(Point3f(100, 50, 100));
    pos.push_back(Point3f(120, 60, 150));
    pos.push_back(Point3f(130, 20, 150));
    face.push_back(0);
    face.push_back(1);
    face.push_back(2);
    face.push_back(3);
    face.push_back(4);
    face.push_back(5);

    //objects.add(CreateTriangleMesh(t, rt, false, 2, 6, pos, face, mat_ptr));
    objects.add(make_shared<Trimesh>(t, rt, false, pos, face, red));*/

    /*std::shared_ptr<Primitive> accelerator =
        CreateBVHAccelerator(std::move(objects.objects));*/

    Transform t = Translate(Point3f(190, 150, 190));
    //cout << t.GetMatrix().m[0][3] << endl;
    Transform rt = Translate(Point3f(-190, -150, -190));

    vector<Point3f> pos;
    vector<int> face;
    char* file = "../ball.obj";
#ifdef USE_OPENMESH
    ReadObjOpenMesh(file, pos, face, 10);
#else
    ReadObj(file, pos, face, 2);
#endif // USE_OPENMESH

    
    objects.add(make_shared<Trimesh>(t, rt, false, pos, face, red));

    /*Transform t1 = Translate(Point3f(400, 150, 190));
    //cout << t.GetMatrix().m[0][3] << endl;
    Transform rt1 = Translate(Point3f(-400, -150, -190));
    objects.add(make_shared<Trimesh>(t1, rt1, false, pos, face, white));*/

    /*std::shared_ptr<Primitive> accelerator =
        CreateKdTreeAccelerator(std::move(objects.objects));*/
    //if (!accelerator) accelerator = std::make_shared<BVHAccel>(primitives);

    /*Transform t = Translate(Point3f(190, 90, 190));
    //cout << t.GetMatrix().m[0][3] << endl;
    Transform rt = Translate(Point3f(-190, -90, -190));

    objects.add(make_shared<sphere>(t,rt,false,90,-90, 90,360, red));*/


    /*std::shared_ptr<Primitive> accelerator =
        CreateBVHAccelerator(std::move(objects.objects));*/
    std::shared_ptr<Primitive> accelerator =
        CreateKdTreeAccelerator(std::move(objects.objects));

    return accelerator;
}

