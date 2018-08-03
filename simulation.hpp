/*
 * SOD model - spore simulation
 *
 * Copyright (C) 2015-2017 by the authors.
 *
 * Authors: Zexi Chen (zchen22 ncsu edu)
 *          Vaclav Petras (wenzeslaus gmail com)
 *          Anna Petrasova (kratochanna gmail com)
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef SPORE_H
#define SPORE_H

#include "raster.hpp"

#include <cmath>
#include <tuple>
#include <random>

// PI is used in the code and M_PI is not guaranteed
// fix it, but prefer the system definition
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#ifndef PI
    #define PI M_PI
#endif

using std::cerr;
using std::endl;

/*! Von Mises Distribution (Circular data distribution)

    mu is the mean angle, expressed in radians between 0 and 2*pi,
    and kappa is the concentration parameter, which must be greater
    than or equal to zero. If kappa is equal to zero, this distribution
    reduces to a uniform random angle over the range 0 to 2*pi.
*/
class von_mises_distribution
{
public:
    von_mises_distribution(double mu, double kappa)
        : mu(mu), kappa(kappa), distribution(0.0, 1.0)
    {}
    template<class Generator>
    double operator ()(Generator& generator)
    {
        double a, b, c, f, r, theta, u1, u2, u3, z;

        if (kappa <= 1.e-06)
            return 2 * PI * distribution(generator);

        a = 1.0 + sqrt(1.0 + 4.0 * kappa * kappa);
        b = (a - sqrt(2.0 * a)) / (2.0 * kappa);
        r = (1.0 + b * b) / (2.0 * b);

        while (true) {
            u1 = distribution(generator);
            z = cos(PI * u1);
            f = (1.0 + r * z) / (r + z);
            c = kappa * (r - f);
            u2 = distribution(generator);
            if (u2 <= c * (2.0 - c) || u2 < c * exp(1.0 - c))
                break;
        }

        u3 = distribution(generator);
        if (u3 > 0.5) {
            theta = fmod(mu + acos(f), 2 * PI);
        }
        else {
            theta = fmod(mu - acos(f), 2 * PI);
        }
        return theta;
    }
private:
    double mu;
    double kappa;
    std::uniform_real_distribution<double> distribution;
};

enum Rtype
{
    CAUCHY, CAUCHY_MIX          // NO means that there is no wind
};

// NONE means that there is no wind
enum Direction
{
    N = 0, NE = 45, E = 90, SE = 135, S = 180, SW = 225, W = 270, NW = 315, NONE
};

/*! The main class to control the spread simulation.
 *
 * The template parameters IntegerRaster and FloatRaster are raster
 * image or matrix types. Any 2D numerical array should work as long as
 * it uses function call operator to access the values, i.e. it provides
 * indexing for reading and writing values using `()`. In other words,
 * the two following operations should be possible:
 *
 * ```
 * a(i, j) = 1;
 * a(i, j) == 1;
 * ```
 *
 * The PoPS library offers a Raster template class to fill this role,
 * but other classes can be used as well.
 */
template<typename IntegerRaster, typename FloatRaster>
class Simulation
{
private:
    int width;
    int height;
    // the west-east resolution of the pixel
    int w_e_res;
    // the north-south resolution of the pixel
    int n_s_res;
    IntegerRaster sp;
    std::default_random_engine generator;
public:

    Simulation(unsigned random_seed, const IntegerRaster &size)
        :
          width(size.getWidth()),
          height(size.getHeight()),
          w_e_res(size.getWEResolution()),
          n_s_res(size.getNSResolution()),
          sp(width, height, w_e_res, n_s_res)
    {
        generator.seed(random_seed);
    }

    Simulation() = delete;

    void SporeRemove(IntegerRaster& I, IntegerRaster& S,
                     const FloatRaster& temperature,
                     double critical_temperature)
    {
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                if (temperature(i, j) < critical_temperature) {
                    S(i, j) += I(i, j);  // move back to suseptable pool
                    I(i, j) = 0;  // remove all infection
                }
            }
        }
    }

    void SporeGen(const IntegerRaster& I, const double *weather,
                  double weather_value, double rate)
    {
        double lambda = 0;
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                if (I(i, j) > 0) {
                    if (weather)
                        lambda = rate * weather[i * width + j];
                    else
                        lambda = rate * weather_value;
                    int sum = 0;
                    std::poisson_distribution<int> distribution(lambda);

                    for (int k = 0; k < I(i, j); k++) {
                        sum += distribution(generator);
                    }
                    sp(i, j) = sum;
                }
                else {
                    sp(i, j) = 0;
                }
            }
        }
    }

    void SporeSpreadDisp_singleSpecies(IntegerRaster& S, IntegerRaster& I, IntegerRaster& I2,
                                       const IntegerRaster& lvtree_rast, std::vector<std::tuple<int, int> > &outside_spores, Rtype rtype,
                                       const double *weather, double weather_value,
                                       double scale1, double kappa = 2,
                                       Direction wdir = NONE, double scale2 = 0.0,
                                       double gamma = 0.0)
    {
        std::cauchy_distribution < double >distribution_cauchy_one(0.0, scale1);
        std::cauchy_distribution < double >distribution_cauchy_two(0.0, scale2);

        std::bernoulli_distribution distribution_bern(gamma);
        std::uniform_real_distribution < double >distribution_uniform(0.0, 1.0);

        if (wdir == NONE)
            kappa = 0;
        von_mises_distribution vonmisesvariate(wdir * PI / 180, kappa);

        double dist = 0;
        double theta = 0;
        bool scale2_used = false;

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                if (sp(i, j) > 0) {
                    for (int k = 0; k < sp(i, j); k++) {

                        // generate the distance from cauchy distribution or cauchy mixture distribution
                        if (rtype == CAUCHY) {
                            dist = abs(distribution_cauchy_one(generator));
                        }
                        else if (rtype == CAUCHY_MIX) {
                            // use bernoulli distribution to act as the sampling with prob(gamma,1-gamma)
                            if (distribution_bern(generator)) {
                                dist = abs(distribution_cauchy_one(generator));
                                scale2_used = false;
                            }
                            else {
                                dist = abs(distribution_cauchy_two(generator));
                                scale2_used = true;
                            }
                        }
                        else {
                            cerr <<
                                    "The paramter Rtype muse be set as either CAUCHY OR CAUCHY_MIX"
                                 << endl;
                            exit(EXIT_FAILURE);
                        }

                        theta = vonmisesvariate(generator);

                        int row = i - round(dist * cos(theta) / n_s_res);
                        int col = j + round(dist * sin(theta) / w_e_res);

                        if (row < 0 || row >= height || col < 0 || col >= width) {
                            // export only spores coming from long-range dispersal kernel outside of modeled area
                            if (scale2_used)
                                outside_spores.emplace_back(std::make_tuple(row, col));
                            continue;
                        }
                        if (S(row, col) > 0) {
                            double prob_S =
                                    (double)(S(row, col)) /
                                    lvtree_rast(row, col);
                            double U = distribution_uniform(generator);

                            if (weather)
                                prob_S *= weather[row * width + col];
                            else
                                prob_S *= weather_value;
                            if (U < prob_S) {
                                I(row, col) += 1;
                                I2(row, col) += 1;
                                S(row, col) -= 1;
                            }
                        }
                    }
                }
            }
        }
    }

    void SporeSpreadDisp(IntegerRaster& S_umca, IntegerRaster& S_oaks, IntegerRaster& I_umca,
                         IntegerRaster& I_oaks, const IntegerRaster& lvtree_rast, Rtype rtype,
                         const double *weather, double weather_value,
                         double scale1, double kappa = 2,
                         Direction wdir = NONE, double scale2 = 0.0,
                         double gamma = 0.0)
    {
        std::cauchy_distribution < double >distribution_cauchy_one(0.0, scale1);
        std::cauchy_distribution < double >distribution_cauchy_two(0.0, scale2);

        std::bernoulli_distribution distribution_bern(gamma);
        std::uniform_real_distribution < double >distribution_uniform(0.0, 1.0);

        if (wdir == NONE)
            kappa = 0;
        von_mises_distribution vonmisesvariate(wdir * PI / 180, kappa);

        double dist = 0;
        double theta = 0;

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                if (sp(i, j) > 0) {
                    for (int k = 0; k < sp(i, j); k++) {

                        // generate the distance from cauchy distribution or cauchy mixture distribution
                        if (rtype == CAUCHY) {
                            dist = abs(distribution_cauchy_one(generator));
                        }
                        else if (rtype == CAUCHY_MIX) {
                            // use bernoulli distribution to act as the sampling with prob(gamma,1-gamma)
                            if (distribution_bern(generator))
                                dist = abs(distribution_cauchy_one(generator));
                            else
                                dist = abs(distribution_cauchy_two(generator));
                        }
                        else {
                            cerr <<
                                    "The paramter Rtype muse be set as either CAUCHY OR CAUCHY_MIX"
                                 << endl;
                            exit(EXIT_FAILURE);
                        }

                        theta = vonmisesvariate(generator);

                        int row = i - round(dist * cos(theta) / n_s_res);
                        int col = j + round(dist * sin(theta) / w_e_res);

                        if (row < 0 || row >= height)
                            continue;
                        if (col < 0 || col >= width)
                            continue;

                        if (row == i && col == j) {
                            if (S_umca(row, col) > 0 ||
                                    S_oaks(row, col) > 0) {
                                double prob =
                                        (double)(S_umca(row, col) +
                                                 S_oaks(row, col)) /
                                        lvtree_rast(row, col);

                                double U = distribution_uniform(generator);

                                if (weather)
                                    prob = prob * weather[row * width + col];
                                else
                                    prob = prob * weather_value;

                                // if U < prob, then one host will become infected
                                if (U < prob) {
                                    double prob_S_umca =
                                            (double)(S_umca(row, col)) /
                                            (S_umca(row, col) +
                                             S_oaks(row, col));
                                    double prob_S_oaks =
                                            (double)(S_oaks(row, col)) /
                                            (S_umca(row, col) +
                                             S_oaks(row, col));

                                    std::bernoulli_distribution
                                        distribution_bern_prob(prob_S_umca);
                                    if (distribution_bern_prob(generator)) {
                                        I_umca(row, col) += 1;
                                        S_umca(row, col) -= 1;
                                    }
                                    else {
                                        I_oaks(row, col) += 1;
                                        S_oaks(row, col) -= 1;
                                    }
                                }
                            }
                        }
                        else {
                            if (S_umca(row, col) > 0) {
                                double prob_S_umca =
                                        (double)(S_umca(row, col)) /
                                        lvtree_rast(row, col);
                                double U = distribution_uniform(generator);

                                if (weather)
                                    prob_S_umca *= weather[row * width + col];
                                else
                                    prob_S_umca *= weather_value;
                                if (U < prob_S_umca) {
                                    I_umca(row, col) += 1;
                                    S_umca(row, col) -= 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

};

// for backwards compatibility
typedef Simulation<Img, DImg> Sporulation;

#endif // SPORE_H
