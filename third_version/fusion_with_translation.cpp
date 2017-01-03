#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <string>

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

//computeWeight(i,j,i+1,j,I1color,I2color,offset1,offset2)
// computes the weight of the edge connecting points (i1,j1) and (i2,j2) in the final image based on their values on coloured images I1 and I2.
//  
double computeWeight(int i1, int j1, int i2, int j2, const Image<Vec3b>&I1color, const Image<Vec3b>&I2color, Point& offset1, Point& offset2){
	Scalar p1I1(I1color(i1-offset1.x,j1-offset1.y));
	Scalar p1I2(I2color(i1-offset2.x,j1-offset2.y));
	Scalar p2I1(I1color(i2-offset1.x,j2-offset1.y));
	Scalar p2I2(I2color(i2-offset2.x,j2-offset2.y));
	return norm(p1I1, p1I2) + norm(p2I1, p2I2);
}

Graph<double,double,double>createGraphFromRectangle(const Rectangle& rec, const Rectangle& overlap, bool right_order1, bool right_order2, const Image<Vec3b>&I1color, const Image<Vec3b>&I2color, Point offset1, Point offset2, int type, int delta){
	Graph<double,double,double> G((rec.p2.x-rec.p1.x)*(rec.p2.y-rec.p1.y),2*(rec.p2.x-rec.p1.x)*(rec.p2.y-rec.p1.y));
	G.add_node((rec.p2.x-rec.p1.x)*(rec.p2.y-rec.p1.y));
	for(int i=rec.p1.x; i<rec.p2.x; i++){
		for(int j=rec.p1.y; j<rec.p2.y; j++){
			// if we are analyzing a point in the intersection region
			if(i>=overlap.p1.x && i<overlap.p2.x && j>=overlap.p1.y && j<overlap.p2.y){
				// we add edges between adjacent points
				if(i<overlap.p2.x-1)
					G.add_edge((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), (i+1-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), computeWeight(i,j,i+1,j,I1color,I2color,offset1,offset2), computeWeight(i,j,i+1,j,I1color,I2color,offset1,offset2));
				// we add edges between adjacent points
				if(j<overlap.p2.y-1)
					G.add_edge((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x), (i-rec.p1.x)+(j+1-rec.p1.y)*(rec.p2.x-rec.p1.x), computeWeight(i,j,i,j+1,I1color,I2color,offset1,offset2), computeWeight(i,j,i,j+1,I1color,I2color,offset1,offset2));

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

void generateImagesFromGraphAndRec(Image<Vec3b>&label, Image<float>&label2, const Graph<double,double,double>&G, const Rectangle& rec, const Rectangle& overlap, bool right_order1, bool right_order2, const Image<Vec3b>&I1color, const Image<Vec3b>&I2color, Point offset1, Point offset2, int type, int delta){
	for (int i=rec.p1.x;i<rec.p2.x;i++)
		for (int j=rec.p1.y;j<rec.p2.y;j++){
			//cout << "i, j = " << i << " " <<  j << " first image dimensions : from " << rec.p1.x << ", " << rec.p1.y << " to " << rec.p2.x << ", " << rec.p2.y << endl;
			label(i-rec.p1.x,j-rec.p1.y)= (G.what_segment((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x)) == Graph<double,double,double>::SOURCE) ? I1color(i-offset1.x,j-offset1.y) : I2color(i-offset2.x,j-offset2.y);
			label2(i-rec.p1.x,j-rec.p1.y)= (G.what_segment((i-rec.p1.x)+(j-rec.p1.y)*(rec.p2.x-rec.p1.x)) == Graph<double,double,double>::SOURCE) ? 1 : 0;
		}
}

Mat mergeTwoImages(const Image<Vec3b>& image1, const Image<Vec3b>& image2){
  int dstWidth = image1.cols;
  int dstHeight = image1.rows * 2;

  Mat dst = Mat(dstHeight, dstWidth, CV_8UC3, cv::Scalar(0,0,0));
  Rect roi(Rect(0,0,image1.cols, image1.rows));
  Mat targetROI = dst(roi);
  image1.copyTo(targetROI);
  targetROI = dst(Rect(0,image1.rows,image1.cols, image1.rows));
  image2.copyTo(targetROI);
  return dst;
}


void do_photomontage(const Image<Vec3b>&I1color, const Image<Vec3b>&I2color, Point offset1, Point offset2, int type=1, int delta=5, bool showCut=false){
	type++;
	bool right_order1=true, right_order2=true;	
	vector<Rectangle>combined_coordinates = rectangleOverlap(I1color, I2color, offset1, offset2, right_order1, right_order2);
	Rectangle overlap, rec;
	selectRectangles(combined_coordinates, rec, overlap, type);
	
	Graph<double,double,double> G = createGraphFromRectangle(rec, overlap, right_order1, right_order2, I1color, I2color, offset1, offset2, type, delta);
	
	double flow=G.maxflow();

	Image<Vec3b> label(rec.p2.x-rec.p1.x,rec.p2.y-rec.p1.y, DataType<Vec3b>::type);
	Image<float> label2(rec.p2.x-rec.p1.x,rec.p2.y-rec.p1.y, DataType<float>::type);
	generateImagesFromGraphAndRec(label, label2, G, rec, overlap, right_order1, right_order2, I1color, I2color, offset1, offset2, type, delta);

	//imshow("Rec"+to_string(type),label);
	//imshow("BlackAndWhite"+to_string(type),label2);
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

int x_1, y_1, x_2, y_2, delta;
int type, showCut;
Image<Vec3b> I1color;
Image<Vec3b> I2color;

void do_pmtg_trackbar(int, void *){
	do_photomontage(I1color, I2color, Point(x_1,y_1), Point(x_2,y_2), type, delta, showCut==1);
}

int main() {
	testGCuts();

	Point offset1(0,85);
	Point offset2(0,0);

	I1color=imread("../img/Raft.jpg");
	I2color=imread("../img/LittleRiver.jpg");
	//imshow("I1",I1color);
	//waitKey();
	//imshow("I2",I2color);
	//waitKey();

	x_1=x_2=y_2=0;
	y_1=85;
	type=1;
	delta=20;
	showCut=0;
	namedWindow("mywindow", WINDOW_AUTOSIZE);
	createTrackbar("Offset x_1", "mywindow", &x_1, 200, do_pmtg_trackbar);
	createTrackbar("Offset y_1", "mywindow", &y_1, 200, do_pmtg_trackbar);
	createTrackbar("Offset x_2", "mywindow", &x_2, 200, do_pmtg_trackbar);
	createTrackbar("Offset y_2", "mywindow", &y_2, 200, do_pmtg_trackbar);
	createTrackbar("Delta", "mywindow", &delta, 200, do_pmtg_trackbar);
	createTrackbar("Type", "mywindow", &type, 1, do_pmtg_trackbar);
	createTrackbar("Image or cut", "mywindow", &showCut, 1, do_pmtg_trackbar);

 
	//do_photomontage(I1color, I2color, offset1, offset2, 1, 20);
	do_photomontage(I1color, I2color, Point(x_1,y_1), Point(x_2,y_2), 1, delta);
	waitKey();

	return 0;
}
