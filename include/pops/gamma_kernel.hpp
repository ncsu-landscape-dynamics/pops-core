/*
 * PoPS model - random uniform dispersal kernel
 *
 * Copyright (C) 2015-2020 by the authors.
 *
 * Authors: Margaret Lawrimore malawrim ncsu edu
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_GAMMA_KERNEL_HPP
#define POPS_GAMMA_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"
#include "lognormal_kernel.hpp"

#include <random>

namespace pops {

using std::pow;
using std::exp;

/*! Dispersal kernel for log secant
 */
class GammaKernel
{
protected:
    double alpha;
    double theta;
    std::gamma_distribution<double> gamma_distribution;

public:
    GammaKernel(double a, double t)
        : alpha(a), theta(t), gamma_distribution(alpha, 1.0 / theta)
    {}

    template<class Generator>
    double random(Generator& generator)
    {
        return std::abs(gamma_distribution(generator));
    }

    double pdf(double x)
    {
        if (x < 0 || alpha < 0 || theta < 0) {
            return 0;
        }
        return 1.0 / (std::tgamma(alpha) * pow(theta, alpha)) * pow(x, (alpha - 1))
               * exp(-x / theta);
    }

    double cdf(double x)
    {
        double sum = 0.0;
        double beta = 1.0 / theta;
        for (int i = 0; i < alpha; i++) {
            // tgamma = (i-1)! used since c++ has no factorial in std lib
            sum += pow(beta * x, i) / std::tgamma(i + 1);
        }
        return 1 - sum * exp(-beta * x);
    }

    // there is no known closed-form solution
    // R and MatLab use an iterative approach (Newton's method)
    double icdf(double x)
    {
        // pick starting approximation using lognormal icdf
        LogNormalKernel lognormal(0, 1);
        double guess = lognormal.icdf(x);
        // TODO add cdf function
        double check = cdf(guess);
        double numiterations = 500;  // will need to adjust this
        double precision = 0.001;  // will need to adjust this
        for (int i = 0; i < numiterations; i++) {
            if (check < (x - precision) || check > (x + precision)) {
                double dif = check - x;
                // if dif is positive guess is greater than needed
                // if dif is negative guess is less than needed
                double past_guess = guess;
                double derivative = (check - x) / pdf(guess);
                // limit size of next guess
                guess = std::max(guess / 10, std::min(guess * 10, guess - derivative));
                check = cdf(guess);
                // Check if we went to far and need to backtrack
                for (int j = 0; j < 10; j++) {
                    if (std::abs(dif) < std::abs(check - x)) {
                        past_guess = guess;
                        guess = (guess + past_guess) / 2.0;
                        check = cdf(guess);
                    }
                }
            }
            else {
                return guess;
            }
        }
        // TODO error message
        std::cout << "unable to find solution to gamma icdf(" << x << ")\n";
        return -1;
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        return type == DispersalKernelType::Gamma;
    }
};

}  // namespace pops

#endif  // POPS_GAMMA_KERNEL_HPP
