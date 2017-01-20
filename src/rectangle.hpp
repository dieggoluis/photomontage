#pragma once

#include <opencv2/imgproc/imgproc.hpp>


typedef struct Rectangle{
	Point p1;
	Point p2;
	Rectangle(int x0, int  y0, int x1, int y1) : p1(x0, y0), p2(x1, y1) {}
} Rectangle;


