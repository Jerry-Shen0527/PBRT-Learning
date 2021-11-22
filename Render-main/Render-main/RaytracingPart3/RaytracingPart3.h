#pragma once

#include <QtWidgets/QWidget>
#include "ui_RaytracingPart3.h"
#include <QPainter>
#include <QImage>
#include <QEvent>
#include "camera.h"
#include "rtweekend.h"
#include "hittable_list.h"
#include "hittable.h"
#include "sphere.h"
#include "material.h"
#include "moving_sphere.h"
#include "aarect.h"
#include "box.h"
//#include "constant_medium.h"
#include "bvh.h"
#include "pdf.h"
#include <iostream>
//#include <omp.h>
#include "spectrum.h"

using namespace std;

class RaytracingPart3 : public QWidget
{
    Q_OBJECT

public:
    RaytracingPart3(QWidget* parent = Q_NULLPTR);
protected:
    void paintEvent(QPaintEvent*);
    void getImage();
   // void getImageChecker();
   // void getImageEmit();

private:
    //color ray_color(const ray& r, const hittable& world, int depth);
    color ray_color(
        const ray& r, const color& background, const hittable& world,
        shared_ptr<hittable>& lights, int depth);
    color ray_color(
        const ray& r, const Spectrum& background, const hittable& world,
        shared_ptr<hittable>& lights, int depth);
    double hit_sphere(const point3& center, double radius, const ray& r);
    vector<int> write_color(color pixel_color, int samples_per_pixel);

    hittable_list cornell_box();
    hittable_list cornell_boxSpe();

    Ui::RaytracingPart3Class ui;
    bool cal_bool = true;
    QString filename;

};
