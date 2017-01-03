#pragma once

#include <vector>
#include <opencv2/imgproc/imgproc.hpp>
#include <iostream>

using namespace cv;
using namespace std;

template <typename T> class Image : public Mat {
public:
	// Constructors
	Image() {}
	Image(const Mat& A):Mat(A) {}
	Image(int w,int h,int type):Mat(h,w,type) {}
	// Accessors
	inline T operator()(int x,int y) const { 
		if(x<0 || x>=cols || y<0 || y>=rows)
			cout << "point " << y << ", " << x << endl;
		return at<T>(y,x); }
	inline T& operator()(int x,int y) { return at<T>(y,x); }
	inline T operator()(const Point& p) const { return at<T>(p.y,p.x); }
	inline T& operator()(const Point& p) { return at<T>(p.y,p.x); }
	//
	inline int width() const { return cols; }
	inline int height() const { return rows; }
	// To display a floating type image
	Image<uchar> greyImage() const {
		double minVal, maxVal;
		minMaxLoc(*this,&minVal,&maxVal);
		Image<uchar> g;
		convertTo(g, CV_8U, 255.0/(maxVal - minVal), -minVal);
		return g;
	}
};

// Harris
vector<Point> harris(const Image<float>& I, double th,int n);
// Correlation
double NCC(const Image<float>& I1,Point m1,const Image<float>& I2,Point m2,int n);
// Correlation with pre-computed means
Image<float> meanImage(const Image<float>& I,int n);
double NCC(const Image<float>& I1,const Image<float>& meanI1,Point m1,const Image<float>& I2,const Image<float>& meanI2,Point m2,int n);


