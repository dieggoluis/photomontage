
#include "rectangleOverlap.h"
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <algorithm>

using namespace std;

#define DEBUG
#ifdef DEBUG
#define debug(x) {cout << #x << " " << x << endl;}
#endif


// returns three rectangles generated from a two translated images

// the first rectangle corresponds to the rectangle with the largest possible width containing points from both images.
// Its leftest points must belong to one image whereas its rightest points must belong to the other image
// in the example, it is the retangle with vertices (0,20), (70, 20), (70,50), (0,50)

// the second rectangle corresponds to the rectangle with the largest possible heigth contaning points from both images.
// Its highest points must belong to one image whereas its rightest points must belong to the other image
// In the example, it is the rectangle with vertices (20, 0), (70, 0), (70, 60), (20, 60)

// the third rectangle corresponds to the intersection of the two rectangles

// position1 and position2 will contain values true if the leftest/toppest points belong to the first rectangle and
// false if the rightest/lowest points belong to the second rectangle
// In the example, it is the rectangle with vertices (20, 20), (70, 20), (70, 50), (20, 50)

/*
        (20,0)
        ________________________
        |                       |
        |                       |
(0, 20) |                       |  
________|_______________________|_____
|       |                       |     |
|       |                       |     |
|       |                       |     |
|       |                       |     |
|       |                       |     |
|       |              (70,50)  |     |
|       |                       |     |
|       |_______________________|     |
|                                     | (90,60)
|_____________________________________|

*/


vector<Rectangle> rectangleOverlap (const Image<Vec3b>& I1, const Image<Vec3b>& I2, Point offset1, Point offset2, bool& position1, bool& position2) {
    pair<Point, Point> i1 (Point(offset1), Point(I1.width()+offset1.x, I1.height()+offset1.y));   
    pair<Point, Point> i2 (Point(offset2), Point(I2.width()+offset2.x, I2.height()+offset2.y));   

    vector<Rectangle> r(3);

    //first possible image
    if(i1.first.x < i2.first.x){
        position1=true;
        r[0].p1 = Point(i1.first.x, max((i1.first).y, (i2.first).y));
        r[0].p2 = Point(i2.second.x, min((i1.second).y, (i2.second).y));
    }
    else if(i1.first.x > i2.first.x){
        position1=false;
        r[0].p1 = Point(i2.first.x, max((i1.first).y, (i2.first).y));
        r[0].p2 = Point(i1.second.x, min((i1.second).y, (i2.second).y));
    }
    else{ // if i1.first.x == i2.first.x
        if(i1.second.x <= i2.second.x){
            position1=true;
            r[0].p1 = Point(i1.first.x, max((i1.first).y, (i2.first).y));
            r[0].p2 = Point(i2.second.x, min((i1.second).y, (i2.second).y));
        }
        else{
            position1=false;
            r[0].p1 = Point(i2.first.x, max((i1.first).y, (i2.first).y));
            r[0].p2 = Point(i1.second.x, min((i1.second).y, (i2.second).y));            
        }
    }
    //r[0].p1 = Point(min((i1.first).x, (i2.first).x), max((i1.first).y, (i2.first).y));
    //r[0].p2 = Point(max((i1.second).x, (i2.second).x), min((i1.second).y, (i2.second).y));

    //second possible image
    if(i1.first.y < i2.first.y){
        position2=true;
        r[1].p1 = Point(max(i1.first.x, i2.first.x), i1.first.y);
        r[1].p2 = Point(min(i1.second.x, i2.second.x), i2.second.y);
    }
    else if(i1.first.y > i2.first.y){
        position2=false;
        r[1].p1 = Point(max(i1.first.x, i2.first.x), i2.first.y);
        r[1].p2 = Point(min(i1.second.x, i2.second.x), i1.second.y);
    }
    else{ // if i1.first.y == i2.first.y
        if(i1.second.y <= i2.second.y){
            position2=true;
            r[1].p1 = Point(max(i1.first.x, i2.first.x), i1.first.y);
            r[1].p2 = Point(min(i1.second.x, i2.second.x), i2.second.y);
        }
        else{
            position2=false;
            r[1].p1 = Point(max(i1.first.x, i2.first.x), i2.first.y);
            r[1].p2 = Point(min(i1.second.x, i2.second.x), i1.second.y);
        }
    }
    //r[1].p1 = Point(max((i1.first).x, (i2.first).x), min((i1.first).y, (i2.first).y));
    //r[1].p2 = Point(min((i1.second).x, (i2.second).x), max((i1.second).y, (i2.second).y));

    //overlap rectangle
    Rectangle r3;
    r[2].p1 = Point(max((i1.first).x, (i2.first).x), max((i1.first).y, (i2.first).y));
    r[2].p2 = Point(min((i1.second).x, (i2.second).x), min((i1.second).y, (i2.second).y));

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

