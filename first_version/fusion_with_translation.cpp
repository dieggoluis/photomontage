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

void do_photomontage(const Image<Vec3b>&I1color, const Image<Vec3b>&I2color, Point offset1, Point offset2, int type=1, int delta=5){
	bool right_order1=true, right_order2=true;
	vector<Rectangle>combined_coordinates = rectangleOverlap(I1color, I2color, offset1, offset2, right_order1, right_order2);
	Rectangle overlap, rec;
	selectRectangles(combined_coordinates, rec, overlap, type);
	
	Graph<double,double,double> G = createGraphFromRectangle(rec, overlap, right_order1, right_order2, I1color, I2color, offset1, offset2, type, delta);
	
	double flow=G.maxflow();

	Image<Vec3b> label(rec.p2.x-rec.p1.x,rec.p2.y-rec.p1.y, DataType<Vec3b>::type);
	Image<float> label2(rec.p2.x-rec.p1.x,rec.p2.y-rec.p1.y, DataType<float>::type);
	generateImagesFromGraphAndRec(label, label2, G, rec, overlap, right_order1, right_order2, I1color, I2color, offset1, offset2, type, delta);

	imshow("Rec"+to_string(type),label);
	imshow("BlackAndWhite"+to_string(type),label2);	
}


int main() {
	testGCuts();

	Image<Vec3b> I1color=imread("../img/Raft.jpg");
	Image<Vec3b> I2color=imread("../img/LittleRiver.jpg");
	imshow("I1",I1color);
	waitKey();
	imshow("I2",I2color);
	waitKey();

	Point offset1(0,85);
	Point offset2(0,0);

	//do_photomontage(I1color, I2color, offset1, offset2, 1, 20);
	do_photomontage(I1color, I2color, offset1, offset2, 2, 20);
	waitKey();

	return 0;
}
