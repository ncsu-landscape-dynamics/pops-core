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

#include <chrono>
#include <random>
#include "Img.h"

#define PI 3.14159265358979323846

enum Rtype
{
    CAUCHY, CAUCHY_MIX          // NO means that there is no wind
};

class Sporulation
{
private:
    double vonmisesvariate(double mu, double kappa);
    int **sp;
    unsigned seed;
    std::default_random_engine generator;
public:
    Sporulation();
    void SporeGen(Img& I, double *weather, double weather_value, double rate);
    void SporeSpreadDisp(Img& S_umca, Img& S_oaks, Img& I_umca,
                         Img& I_oaks, Img& lvtree_rast, Rtype rtype,
                         double *weather, double weather_value,
                         double scale1, double kappa = 2,
                         Direction wdir = NONE, double scale2 = 0.0,
                         double gamma = 0.0);
    ~Sporulation();
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
    int getMonth()
    {
        return month;
    }
    int getYear()
    {
        return year;
    }
    int getDay()
    {
        return day;
    }
};


#endif
