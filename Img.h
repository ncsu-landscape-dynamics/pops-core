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
    int *data;
public:
    Img();
    Img(Img&& other);
    Img(const Img& other);
    //Img(int width,int height);
    Img(const char *fileName);
    Img(int width, int height, int w_e_res, int n_s_res);
    Img& operator=(Img&& other);
    Img& operator=(const Img& other);

    int getWidth() const;
    int getHeight() const;
    int getWEResolution() const;
    int getNSResolution() const;

    const int& operator()(unsigned row, unsigned col) const
    {
        return data[row * width + col];
    }

    int& operator()(unsigned row, unsigned col)
    {
        return data[row * width + col];
    }

    Img operator+(Img & image);
    Img operator-(Img & image);
    Img operator*(int factor);
    ~Img();

    void toGrassRaster(const char *name);
    void toGdal(const char *name, const char *ref_name);

    static Img fromGrassRaster(const char *name);
};

#endif
