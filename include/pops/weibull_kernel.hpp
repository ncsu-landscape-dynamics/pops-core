/*
 * PoPS model - weibull dispersal kernel
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

#ifndef POPS_WEIBULL_KERNEL_HPP
#define POPS_WEIBULL_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"

#include <random>

namespace pops {

using std::pow;
using std::exp;
using std::log;

/*! Dispersal kernel for weibull distribution
 *  class utilized by RadialKernel and DeterministicKernel
 */
class WeibullKernel
{
protected:
    double a;
    double b;
    std::weibull_distribution<double> weibull_distribution;

public:
    WeibullKernel(double a1, double b1) : a(a1), b(b1), weibull_distribution(a, b) {}

    /*!
     *  Returns random value from weibull distribution
     *  Used by RadialKernel to determine location of spread
     *  @param generator uniform random number generator
     *  @return value from weibull distribution
     */
    template<class Generator>
    double random(Generator& generator)
    {
        return std::abs(weibull_distribution(generator));
    }

    /*!
     *  Weibull probability density function
     *  Used by DeterministicKernel to determine location of spread
     *  @param x point within same space of distribution
     *  @return relative likelihood that a random variable would equal x
     */
    double pdf(double x)
    {
        if (x < 0 || b == 0 || a < 0) {
            return 0;
        }
        // Note: If the value inside exp() is too large it returns zero
        // any a >= 2 returns zero
        return ((a / b) * pow(x / b, a - 1) * exp(-pow(x / b, a)));
    }

    /*!
     *  Weibull inverse cumulative distribution (quantile) function
     *  Used by DeterministicKernel to determine maximum distance of spread
     *  @param x proportion of the distribution
     *  @return value in distribution that is less than or equal to probability (x)
     */
    double icdf(double x)
    {
        if (x < 0 || x >= 1 || b <= 0 || a < 0) {
            return 0;
        }
        return a * pow(-(log(1 - x)), (1.0 / b));
    }
};

}  // namespace pops

#endif  // POPS_WEIBULL_KERNEL_HPP
