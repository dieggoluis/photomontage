#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <string>
#include <stdlib.h>

#include "maxflow/graph.h"
#include "image.h"
#include "rectangleOverlap.h"

#define INF numeric_limits<double>::max()/100

using namespace std;

#define DEBUG
#ifdef DEBUG
#define debug(x) {cout << #x << " " << x << endl;}
#endif
// This section shows how to use the library to compute a minimum cut on the following graph:
//
//		        SOURCE
//		       /       \
//		     1/         \6
//		     /      4    \
//		   node0 -----> node1
//		     |   <-----   |
//		     |      3     |
//		     \            /
//		     5\          /1
//		       \        /
//		          SINK
//
///////////////////////////////////////////////////

void testGCuts()
{
	Graph<int,int,int> g(/*estimated # of nodes*/ 2, /*estimated # of edges*/ 1); 
	g.add_node(2); 
	g.add_tweights( 0,   /* capacities */  1, 5 );
	g.add_tweights( 1,   /* capacities */  6, 1 );
	g.add_edge( 0, 1,    /* capacities */  4, 3 );
	int flow = g.maxflow();
	cout << "Flow = " << flow << endl;
	for (int i=0;i<2;i++)
		if (g.what_segment(i) == Graph<int,int,int>::SOURCE)
			cout << i << " is in the SOURCE set" << endl;
		else
			cout << i << " is in the SINK set" << endl;
}
void computeGradient(const Image<Vec3b>& J_0, Image<float>& G, bool blur_image)
{
	int m = J_0.width(), n = J_0.height();

	Mat J;
	if(blur_image)
		GaussianBlur( J_0, J, Size(3,3), 0, 0, BORDER_DEFAULT );
	else
		J = J_0;
	cout << "did bluring" << endl;
	Image<float> I_0(m,n,CV_32F);
	cvtColor(J,I_0,CV_BGR2GRAY);
	cout << "did conversion" << endl;
	int scale = 1;
  int delta = 0;
	int ddepth = CV_16S;
  /// Generate grad_x and grad_y
  Mat grad, grad_x, grad_y;
  Mat abs_grad_x, abs_grad_y;

  /// Gradient X
  //Scharr( I_0, grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
  Sobel(I_0, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
  convertScaleAbs( grad_x, abs_grad_x );
  cout << "gradient x" << endl;
  /// Gradient Y
  //Scharr( I_0, grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
  Sobel( I_0, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );
  convertScaleAbs( grad_y, abs_grad_y );
  cout << "gradient y" << endl;
  /// Total Gradient (approximate)
  addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad );
  cout << "gradient" << endl;
  for(int i=0; i<m; i++)
  	for(int j=0; j<n; j++){
  		//cout << "i,j="<<i<<","<<j<<", m,n="<<m<<","<<n<<endl;
  		if(isnan(grad.at<float>(j,i)))
  			G(i,j)=0;
  		else
  			G(i,j)=grad.at<float>(j,i);
  	}
	cout << "end" << endl;  	

	//cout << "gradx at 1,99: " << abs_grad_x.at<float>(1,99) << endl;
	//cout << "grady at 1,99: " << abs_grad_y.at<float>(1,99) << endl;
  //cout << "gradient at 1,99: " << G(1,99) << endl;

  //imshow("gradient matrix", grad);
  //waitKey();
}



// computeWeight(i,j,i+1,j,lambda,max_lambda,I1color,I2color,G1,G2offset1,offset2)
// computes the weight of the edge connecting points (i1,j1) and (i2,j2) in the final image based on their values on coloured images I1 and I2.
// Consists of a weighted sum of a norm computed using the BGR matrices and a norm computed the matrices of gradients. 
double computeBGRWeight(int i1, int j1, int i2, int j2, const Image<Vec3b>&I1color, const Image<Vec3b>&I2color, Point& offset1, Point& offset2){
	Scalar p1I1(I1color(i1-offset1.x,j1-offset1.y));
	Scalar p1I2(I2color(i1-offset2.x,j1-offset2.y));
	Scalar p2I1(I1color(i2-offset1.x,j2-offset1.y));
	Scalar p2I2(I2color(i2-offset2.x,j2-offset2.y));
	return norm(p1I1, p1I2) + norm(p2I1, p2I2);
}
double computeGradientWeight(int i1, int j1, int i2, int j2, const Image<float>&G1, const Image<float>&G2, Point& offset1, Point& offset2){
	double p1I1 = G1(i1-offset1.x,j1-offset1.y);
	double p1I2 = G2(i1-offset2.x,j1-offset2.y);
	double p2I1 = G1(i2-offset1.x,j2-offset1.y);
	double p2I2 = G2(i2-offset2.x,j2-offset2.y);
	//cout << "i1,j1,i2,j2="<<i1<<","<<j1<<","<<i2<<","<<j2<<endl;
	//cout << "offset1"<<offset1<<",offset2"<<offset2<<endl;
	//cout << p1I1<<","<<p1I2<<","<<p2I1<<","<<p2I2<<endl;
	return abs(p1I1-p1I2) + abs(p2I1-p2I2);	
}
double computeWeight(int i1, int j1, int i2, int j2, int lambda, int max_lambda, const Image<Vec3b>&I1color, const Image<Vec3b>&I2color, const Image<float>&G1, const Image<float>&G2, Point& offset1, Point& offset2){
	double c1 = computeBGRWeight(i1,j1,i2,j2,I1color,I2color,offset1,offset2);
	double c2 = computeGradientWeight(i1,j1,i2,j2,G1,G2,offset1,offset2);
	if(c1>INF || c1<0) c1=INF;
	if(c2>INF || c2<0) c2=INF;
	//cout << "c1="<<c1<<",c2="<<c2<<endl;
	if(max_lambda==0){
		cout << "max_lambda was set to 0, but was supposed to be constant and greater than zero." << endl;
		system("pause");
		return 0;
	}
	double weight;
	if(lambda==0) weight=c1;
	else if(lambda==max_lambda) weight=c2;
	else if(c1>=INF-1||c2>=INF-1)
		weight=INF;
	else{
		weight=( (max_lambda-lambda)*c1 + lambda*c2 )/max_lambda;
		if(weight>=INF-1)
			weight=INF;
	}
	//cout <<"c1,c2="<<c1<<","<<c2<<endl;
	//cout << "weight="<<weight<<endl; 
	//cout << weight << endl;
	return weight;
}

Graph<double,double,double>createGraphFromRectangle(const Rectangle& rec, const Rectangle& overlap, bool right_order1, bool right_order2, const Image<Vec3b>&I1color, const Image<Vec3b>&I2color, const Image<float>&G1, const Image<float>&G2, Point offset1, Point offset2, int type, int delta, int lambda, int max_lambda){
	cout << "lets go" << endl;
	Graph<double,double,double> G((rec.p2.x-rec.p1.x)*(rec.p2.y-rec.p1.y),2*(rec.p2.x-rec.p1.x)*(rec.p2.y-rec.p1.y));
	G.add_node((rec.p2.x-rec.p1.x)*(rec.p2.y-rec.p1.y));
	cout << "go" << endl;
	for(int i=rec.p1.x; i<rec.p2.x; i++){
		for(int j=rec.p1.y; j<rec.p2.y; j++){
			// if we are analyzing a point in the intersection region
			if(i>=overlap.p1.x && i<overlap.p2.x && j>=overlap.p1.y && j<overlap.p2.y){
				// we add edges between adjacent points
				if(i<overlap.p2.x-1)
					G.add_edge((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), (i+1-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), computeWeight(i,j,i+1,j,lambda,max_lambda,I1color,I2color,G1,G2,offset1,offset2), computeWeight(i,j,i+1,j,lambda,max_lambda,I1color,I2color, G1, G2, offset1,offset2));
				// we add edges between adjacent points
				if(j<overlap.p2.y-1)
					G.add_edge((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), (i-rec.p1.x)+(j+1-rec.p1.y)*(rec.p2.x-rec.p1.x), computeWeight(i,j,i,j+1,lambda,max_lambda,I1color,I2color,G1,G2,offset1,offset2), computeWeight(i,j,i,j+1,lambda,max_lambda,I1color,I2color, G1, G2, offset1,offset2));
				// we assign points close to the border to the image they are the closest
				if(type==1 && i<overlap.p1.x+delta){
					if(right_order1)
						G.add_tweights((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), INF, 0);
					else
						G.add_tweights((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), 0, INF);
				}
				if(type==1 && i>overlap.p2.x-delta){
					if(right_order1)
						G.add_tweights((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), 0, INF);
					else
						G.add_tweights((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), INF, 0);
				}
				if(type==2 && j<overlap.p1.y+delta){
					if(right_order2)
						G.add_tweights((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), INF, 0);
					else
						G.add_tweights((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), 0, INF);
				}
				if(type==2 && j>overlap.p2.y-delta){
					if(right_order2)
						G.add_tweights((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), 0, INF);
					else
						G.add_tweights((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), INF, 0);
				}
			}

			// here we assign points that belong to only one of the images to this image
			else if(i<overlap.p1.x || j<overlap.p1.y){
				if(type==1){
					if(right_order1)
						G.add_tweights((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), INF, 0);
					else
						G.add_tweights((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), 0, INF);
				}
				else{
					if(right_order2)
						G.add_tweights((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), INF, 0);
					else
						G.add_tweights((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), 0, INF);
				}
			}
			// here we assign points that belong to only one of the images to this image
			else if(i>=overlap.p2.x || j>=overlap.p2.y){
				if((type==1 && right_order1) || type==2 && right_order2)
						G.add_tweights((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), 0, INF);
				else
					G.add_tweights((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), INF, 0);
			}
			// by logics, we shouldn't enter this case
			else
				cout << "something went wrong" << endl;
		}
	}
	return G;
}
void selectRectangles(const vector<Rectangle>&combined_coordinates, Rectangle& rec, Rectangle& overlap, int type){
	if(type==1)
		rec = combined_coordinates[0]; // coordinates of rectangle in the horizontal
	else if(type==2)
		rec =  combined_coordinates[1]; // coordinates of rectangle in the vertical
	else{
		cout << "wrong type given " << endl;
		return;
	}
	overlap = combined_coordinates[2]; // overlapped rectangle
}

void generateImagesFromGraphAndRec(Image<Vec3b>&label, Image<float>&label2, const Graph<double,double,double>&G, const Rectangle& rec, const Rectangle& overlap, bool right_order1, bool right_order2, const Image<Vec3b>&I1color, const Image<Vec3b>&I2color, Point offset1, Point offset2, int type, int delta, int lambda){
	for (int i=rec.p1.x;i<rec.p2.x;i++)
		for (int j=rec.p1.y;j<rec.p2.y;j++){
			//cout << "i, j = " << i << " " <<  j << " first image dimensions : from " << rec.p1.x << ", " << rec.p1.y << " to " << rec.p2.x << ", " << rec.p2.y << endl;
			label(i-rec.p1.x,j-rec.p1.y)= (G.what_segment((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x)) == Graph<double,double,double>::SOURCE) ? I1color(i-offset1.x,j-offset1.y) : I2color(i-offset2.x,j-offset2.y);
			label2(i-rec.p1.x,j-rec.p1.y)= (G.what_segment((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x)) == Graph<double,double,double>::SOURCE) ? 1 : 0;
		}
}



void do_photomontage(const Image<Vec3b>&I1color, const Image<Vec3b>&I2color, Point offset1, Point offset2, int type=1, int delta=5, bool showCut=false, int lambda=0, int max_lambda=10, bool blur_image=true){
	type++;
	bool right_order1=true, right_order2=true;	
	vector<Rectangle>combined_coordinates = rectangleOverlap(I1color, I2color, offset1, offset2, right_order1, right_order2);
	Rectangle overlap, rec;
	selectRectangles(combined_coordinates, rec, overlap, type);
	cout << "heyhey"<<endl;
	Image<float>G1(I1color.width(), I1color.height(), CV_32F);
	computeGradient(I1color, G1, blur_image);
	cout << "computed first gradient" << endl;
	Image<float>G2(I2color.width(), I2color.height(), CV_32F);
	computeGradient(I2color, G2, blur_image);
	cout << "computed second gradient" << endl;
	Graph<double,double,double> G = createGraphFromRectangle(rec, overlap, right_order1, right_order2, I1color, I2color, G1, G2, offset1, offset2, type, delta, lambda, max_lambda);
	cout << "computed graph" << endl;	
	double flow=G.maxflow();
	cout << "computed flow: " << flow << endl;

	Image<Vec3b> label(rec.p2.x-rec.p1.x,rec.p2.y-rec.p1.y, DataType<Vec3b>::type);
	Image<float> label2(rec.p2.x-rec.p1.x,rec.p2.y-rec.p1.y, DataType<float>::type);
	generateImagesFromGraphAndRec(label, label2, G, rec, overlap, right_order1, right_order2, I1color, I2color, offset1, offset2, type, delta,lambda);
	cout << "generated images" << endl;
	if(showCut){
		imshow("mywindow", label2);
		waitKey();
	}
	else{
		imshow("mywindow", label);
		waitKey();
	}
}
/*
	At first, we use some global variables
	I1color, I2color, x1, y1, x2, y2, delta
*/

int x_1, y_1, x_2, y_2, Delta, Lambda;
const int max_lambda=10;
int Type, ShowCut, Blur_image;
Image<Vec3b> I1color;
Image<Vec3b> I2color;

void do_pmtg_trackbar(int, void *){
	do_photomontage(I1color, I2color, Point(x_1,y_1), Point(x_2,y_2), Type, Delta, ShowCut==1, Lambda, max_lambda, Blur_image);
}

int main() {
	testGCuts();

	Point offset1(0,85);
	Point offset2(0,0);

	I1color=imread("../img/Raft.jpg");
	I2color=imread("../img/LittleRiver.jpg");

	x_1=x_2=y_2=0;
	y_1=85;
	Type=1;
	Delta=20;
	ShowCut=0;
	Lambda=0;
	Blur_image=0;	
	namedWindow("mywindow", WINDOW_AUTOSIZE);
	createTrackbar("Offset x_1", "mywindow", &x_1, 200, do_pmtg_trackbar);
	createTrackbar("Offset y_1", "mywindow", &y_1, 200, do_pmtg_trackbar);
	createTrackbar("Offset x_2", "mywindow", &x_2, 200, do_pmtg_trackbar);
	createTrackbar("Offset y_2", "mywindow", &y_2, 200, do_pmtg_trackbar);
	createTrackbar("Delta", "mywindow", &Delta, 200, do_pmtg_trackbar);
	createTrackbar("Type", "mywindow", &Type, 1, do_pmtg_trackbar);
	createTrackbar("Image/cut", "mywindow", &ShowCut, 1, do_pmtg_trackbar);
	createTrackbar("Pixels/gradient", "mywindow", &Lambda, max_lambda, do_pmtg_trackbar);
	createTrackbar("Blur image", "mywindow", &Blur_image, 1, do_pmtg_trackbar);
 
	//do_photomontage(I1color, I2color, offset1, offset2, 1, 20);
	do_photomontage(I1color, I2color, Point(x_1,y_1), Point(x_2,y_2), 1, Delta,false,Lambda,max_lambda,true);
	waitKey();

	return 0;
}
