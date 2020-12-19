#ifdef POPS_TEST

/*
 * Tests if random number generators appropriately approximates the power law
 * distribution and the hyperbolic secant distribution
 *
 * Copyright (C) 2018-2020 by the authors.
 *
 * Authors: Margaret Lawrimore <malawrim ncsu edu>
 *
 * This file is part of PoPS.

 * PoPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * PoPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with PoPS. If not, see <https://www.gnu.org/licenses/>.
 */

#include <pops/power_law_kernel.hpp>
#include <pops/hyperbolic_secant_kernel.hpp>
#include <pops/logistic_kernel.hpp>
#include <pops/exponential_power_kernel.hpp>
#include <pops/gamma_kernel.hpp>
#include <pops/deterministic_kernel.hpp>
#include <iostream>
#include <fstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef PI
#define PI M_PI
#endif

using std::string;
using std::cout;

using namespace pops;

// These tests are just being used to create a .csv of values produced by the
// RNG that I will then plot and test against a chi square goodness of
// fit test to make sure it aligns appropriately with the desired distribution
int test_power_law_rng()
{
    // TODO xmin has to be the lower bound for x values
    // so if range has to be 0 - 1 then xmin always needs to be 0
    std::default_random_engine generator;
    // Pulling code out of power_law_kernel to test here to be able to save the
    // value produced by the rng
    std::ofstream file;
    file.open("testing_power_law.csv");
    file << "x,y,\n";

    double alpha = 1.5;
    double xmin = 0.01;
    PowerLawKernel power_law_distribution(alpha, xmin);
    std::uniform_real_distribution<double> distribution(0.1, 1.0);
    for (int i = 0; i < 10000; i++) {
        double x = distribution(generator);
        double y = power_law_distribution.icdf(x);
        file << x << "," << y << ",\n";
    }

    // if you expect x% to be between 1-2 then compare that with the actual number
    // TODO number of buckets?
    file.close();
    return 0;
}
int test_hyperbolic_secant_rng()
{
    std::default_random_engine generator;
    double sigma_ = 2.0;
    // HyperbolicSecantDistributionRandom logsech(1.0, 2.0);
    // std::cout << "Hyperbolic secant random: " << logsech.random(generator) << "\n";
    // Pulling code out of hyperbolic_secant_kernel to test here to be able to save the
    // value produced by the rng
    HyperbolicSecantKernel hyperbolic_secant_distribution(sigma_);
    std::ofstream file;
    file.open("testing_secant.csv");
    file << "x,y,\n";

    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (int i = 0; i < 10000; i++) {
        double x = distribution(generator);
        // get random value from a uniform distribution and use it
        // to get a random value from the distribution
        double y = hyperbolic_secant_distribution.icdf(x);
        file << x << "," << y << ",\n";
    }
    file.close();
    return 0;
}
int test_logistic_rng()
{
    std::default_random_engine generator;
    double s_ = 2.0;
    LogisticKernel logistic_distribution(s_);
    std::ofstream file;
    file.open("testing_logistic.csv");
    file << "x,y,\n";

    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (int i = 0; i < 10000; i++) {
        double x = distribution(generator);

        // get random value from a uniform distribution and use it
        // to get a random value from the distribution
        double y = logistic_distribution.icdf(x);
        file << x << "," << y << ",\n";
    }
    file.close();
    return 0;
}

int test_exponential_power_rng()
{
    std::default_random_engine generator;
    double alpha_ = 1.5;
    double beta_ = 0.5;
    ExponentialPowerKernel exponential_power_distribution(alpha_, beta_);
    std::ofstream file;
    file.open("testing_exponential_power.csv");
    file << "x,y,\n";

    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (int i = 0; i < 10000; i++) {
        double x = distribution(generator);
        double y = exponential_power_distribution.icdf(x);
        file << x << "," << y << ",\n";
    }
    file.close();
    return 0;
}

int test_gamma()
{
    std::default_random_engine generator;
    double alpha_ = 1.5;
    double theta_ = 0.5;
    GammaKernel gamma_distribution(alpha_, theta_);
    std::ofstream file;
    file.open("testing_gamma.csv");
    file << "x,y,\n";
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (int i = 0; i < 10000; i++) {
        double x = distribution(generator);
        double y = gamma_distribution.icdf(x);
        file << x << "," << y << ",\n";
    }
    file.close();
    return 0;
}
int main()
{
    int ret = 0;

    ret += test_power_law_rng();
    ret += test_hyperbolic_secant_rng();
    ret += test_logistic_rng();
    ret += test_exponential_power_rng();
    ret += test_gamma();

    std::cout << "Test deterministic number of errors: " << ret << std::endl;
    return ret;
}
#endif  // POPS_TEST
