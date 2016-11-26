
#include "rectangleOverlap.h"
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <algorithm>

using namespace std;

#define DEBUG
#ifdef DEBUG
#define debug(x) {cout << #x << " " << x << endl;}
#endif

vector<Rectangle> rectangleOverlap (Image<Vec3b>& I1, Image<Vec3b>& I2, Point offset1, Point offset2, bool& position1, bool& position2) {
    pair<Point, Point> i1 (Point(offset1), Point(I1.width()+offset1.x, I1.height()+offset1.y));   
    pair<Point, Point> i2 (Point(offset2), Point(I2.width()+offset2.x, I2.height()+offset2.y));   

    vector<Rectangle> r(3);

    //first possible image
    r[0].p1 = Point(min((i1.first).x, (i2.first).x), max((i1.first).y, (i2.first).y));
    r[0].p2 = Point(max((i1.second).x, (i2.second).x), min((i1.second).y, (i2.second).y));

    //second possible image
    Rectangle r2;
    r[1].p1 = Point(max((i1.first).x, (i2.first).x), min((i1.first).y, (i2.first).y));
    r[1].p2 = Point(min((i1.second).x, (i2.second).x), max((i1.second).y, (i2.second).y));

    //overlap triangle
    Rectangle r3;
    r[2].p1 = Point(max((i1.first).x, (i2.first).x), max((i1.first).y, (i2.first).y));
    r[2].p2 = Point(min((i1.second).x, (i2.second).x), min((i1.second).y, (i2.second).y));

    position1 = (r[0].p1).x == (i1.first).x;
    position2 = (r[0].p1).y == (i1.first).y;

    return r;
}

//#define TEST
#ifdef TEST

int main (int argc, char *argv[]) {
    Image<Vec3b> Icolor1 = imread(argv[1]);
    Image<Vec3b> Icolor2 = imread(argv[2]);
   
    imshow("image1", Icolor1);
    imshow("image2", Icolor2);
    debug(Icolor1.cols);
    debug(Icolor1.rows);

    debug(Icolor2.cols);
    debug(Icolor2.rows);
    
    bool pos1, pos2;
    vector<Rectangle> r =  rectangleOverlap(Icolor1, Icolor2, Point(0,0), Point(0,0), pos1, pos2);
    
    debug(r[0].p1.x);
    debug(r[0].p1.y);
    debug(r[0].p2.x);
    debug(r[0].p2.y);

    debug(r[1].p1.x);
    debug(r[1].p1.y);
    debug(r[1].p2.x);
    debug(r[1].p2.y);

    debug(r[2].p1.x);
    debug(r[2].p1.y);
    debug(r[2].p2.x);
    debug(r[2].p2.y);

    waitKey();
    return 0;
}

#endif

