#pragma once

#include "image.h"
#include <opencv2/imgproc/imgproc.hpp>

typedef struct Rectangle {
    Point p1;
    Point p2;   
} Rectangle;

vector<Rectangle> rectangleOverlap (const Image<Vec3b>& I1, const Image<Vec3b>& I2, Point offset1, Point offset2, bool& position1, bool& position2);
