#include "image.h"

// Harris points
vector<Point> harris(const Image<float>& I, double th,int n) {
	vector<Point> v;
	Image<float> H;
	cornerHarris(I,H,10,3,0.04);
	for (int y=n;y<I.height()-n;y++)
		for (int x=n;x<I.width()-n;x++)
			if (H(x,y) > th
				&& H(x,y)>H(x,y+1) && H(x,y)>H(x,y-1) && H(x,y)>H(x-1,y-1) && H(x,y)>H(x-1,y)
				&& H(x,y)>H(x-1,y+1) && H(x,y)>H(x+1,y-1) && H(x,y)>H(x+1,y) && H(x,y)>H(x+1,y+1))
				v.push_back(Point(x,y));
	return v;
}

// Correlation
double mean(const Image<float>& I,Point m,int n) {
	double s=0;
	for (int j=-n;j<=n;j++)
		for (int i=-n;i<=n;i++) 
			s+=I(m+Point(i,j));
	return s/(2*n+1)/(2*n+1);
}

double corr(const Image<float>& I1,Point m1,const Image<float>& I2,Point m2,int n) {
	double M1=mean(I1,m1,n);
	double M2=mean(I2,m2,n);
	double rho=0;
	for (int j=-n;j<=n;j++)
		for (int i=-n;i<=n;i++) {
			rho+=(I1(m1+Point(i,j))-M1)*(I2(m2+Point(i,j))-M2);
		}
		return rho;
}

double NCC(const Image<float>& I1,Point m1,const Image<float>& I2,Point m2,int n) {
	if (m1.x<n || m1.x>=I1.width()-n || m1.y<n || m1.y>=I1.height()-n) return -1;
	if (m2.x<n || m2.x>=I2.width()-n || m2.y<n || m2.y>=I2.height()-n) return -1;
	double c1=corr(I1,m1,I1,m1,n);
	if (c1==0) return -1;
	double c2=corr(I2,m2,I2,m2,n);
	if (c2==0) return -1;
	return corr(I1,m1,I2,m2,n)/sqrt(c1*c2);
}

// ===========================================================

// Correlation with pre-computed means
Image<float> meanImage(const Image<float>& I,int n) {
	Image<float> meanI(I.width(),I.height(),CV_32F);
	for (int j=n;j<I.height()-n;j++) {
		for (int i=n;i<I.width()-n;i++) {
			double s=0;
			for (int dj=-n;dj<=n;dj++)
				for (int di=-n;di<=n;di++) 
					s+=I(i+di,j+dj);
			meanI(i,j)=float(s/(2*n+1)/(2*n+1));
		}
	}
	return meanI;
}

double corr(const Image<float>& I1,const Image<float>& meanI1,Point m1,const Image<float>& I2,const Image<float>& meanI2,Point m2,int n) {
	double M1=meanI1(m1);
	double M2=meanI2(m2);

	double rho=0;
	for (int j=-n;j<=n;j++)
		for (int i=-n;i<=n;i++) {
			rho+=(I1(m1+Point(i,j))-M1)*(I2(m2+Point(i,j))-M2);
		}
	return rho;
}

double NCC(const Image<float>& I1,const Image<float>& meanI1,Point m1,const Image<float>& I2,const Image<float>& meanI2,Point m2,int n) {
	if (m1.x<n || m1.x>=I1.width()-n || m1.y<n || m1.y>=I1.height()-n) return -1;
	if (m2.x<n || m2.x>=I2.width()-n || m2.y<n || m2.y>=I2.height()-n) return -1;
	double c1=corr(I1,meanI1,m1,I1,meanI1,m1,n);
	if (c1==0) return -1;
	double c2=corr(I2,meanI2,m2,I2,meanI2,m2,n);
	if (c2==0) return -1;
	return corr(I1,meanI1,m1,I2,meanI2,m2,n)/sqrt(c1*c2);
}

