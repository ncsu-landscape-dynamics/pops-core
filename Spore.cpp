/*
 *
 *  SOD-model-cpp
 *
 *  Created on: Oct,2015
 *  Author: Zexi Chen(zchen22@ncsu.edu)
 *
 *
 */

#include "Spore.h"

#include <cmath>

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


Sporulation::Sporulation(unsigned random_seed)
{
    this->sp = NULL;
    //seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(random_seed);
}

void Sporulation::SporeGen(Img & I, double *weather, double weather_value, double rate)
{

    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //std::default_random_engine generator(seed);

    int height = I.getHeight();
    int width = I.getWidth();

    if (!sp) {
        sp = (int **)std::malloc(sizeof(int *) * height);
        int *stream = (int *)std::malloc(sizeof(int) * width * height);


        for (int i = 0; i < height; i++) {
            sp[i] = &stream[i * width];
        }
    }

    double lambda = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (I.data[i][j] > 0) {
                if (weather)
                    lambda = rate * weather[i * width + j];
                else
                    lambda = rate * weather_value;
                int sum = 0;
                std::poisson_distribution<int> distribution(lambda);

                for (int k = 0; k < I.data[i][j]; k++) {
                    sum += distribution(generator);
                }
                sp[i][j] = sum;
            }
            else {
                sp[i][j] = 0;
            }
        }
    }
}

void Sporulation::SporeSpreadDisp(Img & S_umca, Img & S_oaks, Img & I_umca,
                                  Img & I_oaks, Img & lvtree_rast,
                                  Rtype rtype, double *weather,
                                  double weather_value, double scale1,
                                  double kappa, Direction wdir, double scale2,
                                  double gamma)
{

    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //std::default_random_engine generator(seed);
    std::cauchy_distribution < double >distribution_cauchy_one(0.0, scale1);
    std::cauchy_distribution < double >distribution_cauchy_two(0.0, scale2);

    std::bernoulli_distribution distribution_bern(gamma);
    std::uniform_real_distribution < double >distribution_uniform(0.0, 1.0);


    int height = S_umca.getHeight();
    int width = S_umca.getWidth();
    int w_e_res = S_umca.getWEResolution();
    int n_s_res = S_umca.getNSResolution();

    double dist = 0;
    double theta = 0;

    if (!sp) {
        cerr << "The spore matrix is empty!" << endl;
        return;
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (sp[i][j] > 0) {
                for (int k = 0; k < sp[i][j]; k++) {

                    // generate the distance from cauchy distribution or cauchy mixture distribution
                    if (rtype == CAUCHY) {
                        dist = abs(distribution_cauchy_one(generator));
                    }
                    else if (rtype == CAUCHY_MIX) {
                        if (gamma >= 1 || gamma <= 0) {
                            cerr <<
                                    "The parameter gamma must be in the range (0~1)"
                                 << endl;
                            return;
                        }
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

                    if (wdir == NONE) {
                        kappa = 0;
                    }

                    theta = vonmisesvariate(wdir * PI / 180, kappa);

                    int row = i - round(dist * cos(theta) / n_s_res);
                    int col = j + round(dist * sin(theta) / w_e_res);

                    if (row < 0 || row >= height)
                        continue;
                    if (col < 0 || col >= width)
                        continue;

                    if (row == i && col == j) {
                        if (S_umca.data[row][col] > 0 ||
                                S_oaks.data[row][col] > 0) {
                            double prob =
                                    (double)(S_umca.data[row][col] +
                                             S_oaks.data[row][col]) /
                                    lvtree_rast.data[row][col];

                            double U = distribution_uniform(generator);

                            if (weather)
                                prob = prob * weather[row * width + col];
                            else
                                prob = prob * weather_value;

                            // if U < prob, then one host will become infected
                            if (U < prob) {
                                double prob_S_umca =
                                        (double)(S_umca.data[row][col]) /
                                        (S_umca.data[row][col] +
                                         S_oaks.data[row][col]);
                                double prob_S_oaks =
                                        (double)(S_oaks.data[row][col]) /
                                        (S_umca.data[row][col] +
                                         S_oaks.data[row][col]);

                                std::bernoulli_distribution
                                    distribution_bern_prob(prob_S_umca);
                                if (distribution_bern_prob(generator)) {
                                    I_umca.data[row][col] += 1;
                                    S_umca.data[row][col] -= 1;
                                }
                                else {
                                    I_oaks.data[row][col] += 1;
                                    S_oaks.data[row][col] -= 1;
                                }
                            }
                        }
                    }
                    else {
                        if (S_umca.data[row][col] > 0) {
                            double prob_S_umca =
                                    (double)(S_umca.data[row][col]) /
                                    lvtree_rast.data[row][col];
                            double U = distribution_uniform(generator);

                            if (weather)
                                prob_S_umca *= weather[row * width + col];
                            else
                                prob_S_umca *= weather_value;
                            if (U < prob_S_umca) {
                                I_umca.data[row][col] += 1;
                                S_umca.data[row][col] -= 1;
                            }
                        }
                    }
                }
            }
        }
    }
}

double Sporulation::vonmisesvariate(double mu, double kappa)
{

    /**
    Von Mises Distribution(Circular data distribution)

    mu is the mean angle, expressed in radians between 0 and 2*pi,
    and kappa is the concentration parameter, which must be greater
    than or equal to zero. If kappa is equal to zero, this distribution
    reduces to a uniform random angle over the range 0 to 2*pi

    */

    double a, b, c, f, r, theta, u1, u2, u3, z;
    std::uniform_real_distribution < double >distribution(0.0, 1.0);

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

Sporulation::~Sporulation()
{
    if (sp) {
        if (sp[0])
            delete[]sp[0];
        delete[]sp;
    }
}

bool Date::compareDate(Date & endtime)
{
    if (year < endtime.year) {
        return true;
    }
    else if (year == endtime.year) {
        if (month < endtime.month) {
            return true;
        }
        else if (month == endtime.month) {
            if (day <= endtime.day)
                return true;
            else
                return false;
        }
        else {
            return false;
        }
    }
    else {
        return false;
    }
}

void Date::increasedByWeek()
{
    day += 7;
    if (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0)) {
        if (day > day_in_month[1][month]) {
            day = day - day_in_month[1][month];
            month++;
            if (month > 12) {
                year++;
                month = 1;
            }
        }
    }
    else {
        if (day > day_in_month[0][month]) {
            day = day - day_in_month[0][month];
            month++;
            if (month > 12) {
                year++;
                month = 1;
            }
        }
    }
}
