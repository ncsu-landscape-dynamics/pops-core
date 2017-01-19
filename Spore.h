/*
 *
 *  SOD-model-cpp
 *
 *  Created on: Oct,2015
 *  Author: Zexi Chen(zchen22@ncsu.edu)
 *
 *
 */


#ifndef SPORE_H
#define SPORE_H

#include "Img.h"

#include <chrono>
#include <random>


enum Rtype
{
    CAUCHY, CAUCHY_MIX          // NO means that there is no wind
};

class Sporulation
{
private:
    double vonmisesvariate(double mu, double kappa);
    int width;
    int height;
    // the west-east resolution of the pixel
    int w_e_res;
    // the north-south resolution of the pixel
    int n_s_res;
    Img sp;
    std::default_random_engine generator;
public:
    Sporulation(unsigned random_seed, const Img &size);
    void SporeGen(Img& I, double *weather, double weather_value, double rate);
    void SporeSpreadDisp(Img& S_umca, Img& S_oaks, Img& I_umca,
                         Img& I_oaks, Img& lvtree_rast, Rtype rtype,
                         double *weather, double weather_value,
                         double scale1, double kappa = 2,
                         Direction wdir = NONE, double scale2 = 0.0,
                         double gamma = 0.0);
};

class Date
{
private:
    int year;
    int month;
    int day;
    int day_in_month[2][13] = {
        {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
        {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
    };

public:
    Date(int y, int m, int d):year(y), month(m), day(d)
    {
    }
    bool compareDate(Date & endtime);
    void increasedByWeek();
    int getMonth() const
    {
        return month;
    }
    int getYear() const
    {
        return year;
    }
    int getDay() const
    {
        return day;
    }
};


#endif
