/*
 * PoPS model - hyperbolic secant dispersal kernel
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

#ifndef POPS_HYPERBOLIC_SECANT_KERNEL_HPP
#define POPS_HYPERBOLIC_SECANT_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"

#include <random>

namespace pops {

using std::cosh;
using std::tan;
using std::log;

/*! Dispersal kernel for hyperbolic secant
 *  class utilized by RadialKernel and DeterministicKernel
 */
class HyperbolicSecantKernel
{
protected:
    double sigma;

public:
    HyperbolicSecantKernel(double s) : sigma(s) {}

    /*!
     *  Returns random value from hyperbolic secant distribution
     *  Used by RadialKernel to determine location of spread
     *  @param generator uniform random number generator
     *  @return value from hyperbolic secant distribution
     */
    template<class Generator>
    double random(Generator& generator)
    {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        double x = distribution(generator);
        // get random value from a uniform distribution and use it
        // to get a random value from the distribution
        // if (mu_ != 0 && sigma_ != 1) {
        //    return ((log(tan((x * M_PI) / 2.0)) * (2.0 * sigma_)) / M_PI) + mu_;
        //}
        // return (2 / M_PI) * log(tan(M_PI / 2 * x));
        return icdf(x);
    }

    /*!
     *  Hyperbolic secant probability density function
     *  Used by DeterministicKernel to determine location of spread
     *  @param x point within same space of distribution
     *  @return relative likelihood that a random variable would equal x
     */
    double pdf(double x)
    {
        if (x <= 0 || sigma == 0) {
            return 0;
        }
        if (sigma == 1) {
            return 0.5 * (1 / cosh((M_PI * x) / 2));
        }
        return (1.0 / (2 * sigma)) * (1 / cosh((M_PI * x) / (2 * sigma)));
    }

    /*!
     *  Hyperbolic secant inverse cumulative distribution (quantile) function
     *  Used by DeterministicKernel to determine maximum distance of spread
     *  @param x proportion of the distribution
     *  @return value in distribution that is less than or equal to probability (x)
     */
    double icdf(double x)
    {
        if (x <= 0 || sigma == 0) {
            return 0;
        }
        if (sigma == 1) {
            return (2.0 / M_PI) * log(tan(M_PI / 2.0 * x));
        }
        return ((log(tan((x * M_PI) / 2.0)) * (2.0 * sigma)) / M_PI);
    }
};

}  // namespace pops

#endif  // POPS_HYPERBOLIC_SECANT_KERNEL_HPP
