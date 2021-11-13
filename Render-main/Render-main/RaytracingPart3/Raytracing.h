#pragma once

#include <QtWidgets/QWidget>
#include "ui_Raytracing.h"
#include <QPainter>
#include <QImage>
#include <QEvent>
#include "camera.h"
#include "rtweekend.h"
#include "hittable_list.h"
#include "sphere.h"
#include "material.h"
#include "moving_sphere.h"
#include "aarect.h"
#include "box.h"
#include "constant_medium.h"
#include "bvh.h"

using namespace std;

class Raytracing : public QWidget
{
    Q_OBJECT

public:
    Raytracing(QWidget *parent = Q_NULLPTR);
protected:
    void paintEvent(QPaintEvent*);
    void getImage();
    void getImageChecker();
    void getImageEmit();
    
private:
    color ray_color(const ray& r, const hittable& world, int depth);
    color ray_color(const ray& r, const color& background, const hittable& world, int depth);
    double hit_sphere(const point3& center, double radius, const ray& r);
    vector<int> write_color(color pixel_color, int samples_per_pixel);

    hittable_list random_scene();
    hittable_list two_spheres();
    hittable_list two_perlin_spheres();
    hittable_list earth();
    hittable_list simple_light();
    hittable_list cornell_box();
    hittable_list cornell_smoke();
    hittable_list final_scene();

    Ui::RaytracingClass ui;
    bool cal_bool = true;
    QString filename;
    
};
