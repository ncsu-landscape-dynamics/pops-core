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
    void SporeGen(const Img& I, const double *weather,
                  double weather_value, double rate);
    void SporeSpreadDisp(Img& S_umca, Img& S_oaks, Img& I_umca,
                         Img& I_oaks, const Img& lvtree_rast, Rtype rtype,
                         const double *weather, double weather_value,
                         double scale1, double kappa = 2,
                         Direction wdir = NONE, double scale2 = 0.0,
                         double gamma = 0.0);
};

#endif
