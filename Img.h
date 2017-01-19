/*
 *
 *  SOD-model-cpp
 *
 *  Created on: Oct,2015
 *  Author: Zexi Chen(zchen22@ncsu.edu)
 *
 *
 */


#ifndef IMG_H
#define IMG_H

#include <iostream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <ctime>


enum Direction
{
    N = 0, NE = 45, E = 90, SE = 135, S = 180, SW = 225, W = 270, NW = 315, NONE  // NO means that there is no wind
};

class Img
{
private:
    int width;
    int height;
    // the west-east resolution of the pixel
    int w_e_res;
    // the north-south resolution of the pixel
    int n_s_res;
public:
    int **data;
    Img();
    //Img(int width,int height);
    Img(const char *fileName);
    Img(int width, int height, int w_e_res, int n_s_res, int **data);
    int getWidth();
    int getHeight();
    int getWEResolution();
    int getNSResolution();

    Img operator+(Img & image);
    Img operator-(Img & image);
    Img operator*(int factor);
    ~Img();

    void toGrassRaster(const char *name);
    void toGdal(const char *name, const char *ref_name);

    static Img fromGrassRaster(const char *name);
};

#endif
