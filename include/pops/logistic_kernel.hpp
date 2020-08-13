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

#ifndef POPS_LOGISTIC_KERNEL_HPP
#define POPS_LOGISTIC_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"

#include <random>

namespace pops {

using std::pow;

/*! Dispersal kernel for log secant
 */
class LogisticKernel
{
protected:
    double mu;
    double s;

public:
    LogisticKernel(double m, double ss) : mu(m), s(ss) {}

    template<class Generator>
    double random(Generator& generator)
    {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        double x = distribution(generator);
        // get random value from a uniform distribution and use it
        // to get a random value from the distribution
        // if (mu_ == 0 && s_ == 1) {
        //    return log(x / (1 - x));
        //}
        // return mu_ + s_ * log(x / (1 - x));
        return icdf(x);
    }

    double pdf(double x)
    {
        if (x < 0 || s == 0) {
            return 0;
        }
        if (mu == 0 && s == 1) {
            return exp(-x) / pow(1 + exp(-x), 2);
        }
        return (exp(-(x - mu) / s)) / (s * pow(1 + exp(-(x - mu) / s), 2));
    }

    double icdf(double x)
    {
        if (x <= 0 || s == 0) {
            return 0;
        }
        if (mu == 0 && s == 1) {
            return log(x / (1.0 - x));
        }
        return mu + (s * log(x / (1.0 - x)));
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        return type == DispersalKernelType::Logistic;
    }
};

}  // namespace pops

#endif  // POPS_LOGISTIC_KERNEL_HPP