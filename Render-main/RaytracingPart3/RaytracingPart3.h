#pragma once

#include <QtWidgets/QWidget>
#include "ui_RaytracingPart3.h"
#include <QPainter>
#include <QImage>
#include <QEvent>
#include "camera.h"
#include "rtweekend.h"
#include "transform.h"
#include "primitive.h"
#include "primitive_list.h"
#include "sphere.h"
#include "material.h"
//#include "moving_sphere.h"
#include "aarect.h"
//#include "box.h"
//#include "constant_medium.h"
//#include "bvh.h"
#include "trianglenew.h"
#include "pdf.h"
#include <iostream>
//#include <omp.h>
#include "spectrum.h"
#include "kdtreeaccel.h"
#include "bvh.h"
#include "trimesh.h"

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
    //color ray_color(const ray& r, const Primitive& world, int depth);
    color ray_color(
        const ray& r, const color& background, const Primitive& world,
        shared_ptr<Primitive>& lights, int depth);
    color ray_color(
        const ray& r, const color& background, const Primitive& world,
        int depth);
    color ray_color(
        const ray& r, const Spectrum& background, const Primitive& world,
        shared_ptr<Primitive>& lights, int depth);
    double hit_sphere(const Point3f& center, double radius, const ray& r);
    vector<int> write_color(color pixel_color, int samples_per_pixel);

    //Primitive_list cornell_box();
    //Primitive_list cornell_boxSpe();
    Primitive_list Test();
    

    Ui::RaytracingPart3Class ui;
    bool cal_bool = true;
    QString filename;

};
