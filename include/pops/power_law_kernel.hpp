/*
 * PoPS model - power law dispersal kernel
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

#ifndef POPS_POWER_LAW_KERNEL_HPP
#define POPS_POWER_LAW_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"

#include <random>

namespace pops {

using std::pow;

/*! Dispersal kernel for power law distribution
 *  class utilized by RadialKernel and DeterministicKernel
 */
class PowerLawKernel
{
protected:
    double xmin;
    double alpha;

public:
    PowerLawKernel(double xm, double a) : xmin(xm), alpha(a) {}

    /*!
     *  Returns random value from power law distribution
     *  Used by RadialKernel to determine location of spread
     *  @param generator uniform random number generator
     *  @return value from power law distribution
     */
    template<class Generator>
    double random(Generator& generator)
    {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        double x = distribution(generator);
        // since distribution is 0 to 1 xmin needs to be greater than 0
        if (xmin <= 0) {
            return 0;
        }
        return icdf(x);
    }

    /*!
     *  Power law probability density function
     *  Used by DeterministicKernel to determine location of spread
     *  @param x point within same space of distribution
     *  @return relative likelihood that a random variable would equal x
     *  @note only works with alpha < 1 so center square in matrix is always zero
     */
    double pdf(double x)
    {
        if (x <= 0 || xmin == 0 || alpha <= 1.0) {
            return 0;
        }
        return ((alpha - 1.0) / xmin) * pow(x / xmin, -alpha);
    }

    /*!
     *  Power law inverse cumulative distribution (quantile) function
     *  Used by DeterministicKernel to determine maximum distance of spread
     *  @param x proportion of the distribution
     *  @return value in distribution that is less than or equal to probability (x)
     */
    double icdf(double x)
    {
        if (x <= 0 || x >= 1 || xmin == 0 || alpha <= 1.0) {
            return 0;
        }
        return pow(x, (1.0 / (-alpha + 1.0))) * xmin;
    }
};

}  // namespace pops

#endif  // POPS_POWER_LAW_KERNEL_HPP
