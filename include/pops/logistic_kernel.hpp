/*
 * PoPS model - logistic dispersal kernel
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
using std::exp;
using std::log;
/*! Dispersal kernel for logistic distribution
 *  class utilized by RadialKernel and DeterministicKernel
 */
class LogisticKernel
{
protected:
    double s;

public:
    LogisticKernel(double scale, double unused) : s(scale) {}

    /*!
     *  Returns random value from logistic distribution
     *  Used by RadialKernel to determine location of spread
     *  @param generator uniform random number generator
     *  @return value from logistic distribution
     */
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

    /*!
     *  Logistic probability density function
     *  Used by DeterministicKernel to determine location of spread
     *  @param x point within same space of distribution
     *  @return relative likelihood that a random variable would equal x
     */
    double pdf(double x)
    {
        if (x < 0 || s == 0) {
            return 0;
        }
        if (s == 1) {
            return exp(-x) / pow(1 + exp(-x), 2);
        }
        return (exp(-x / s)) / (s * pow(1 + exp(-x / s), 2));
    }

    /*!
     *  Logistic inverse cumulative distribution (quantile) function
     *  Used by DeterministicKernel to determine maximum distance of spread
     *  @param x proportion of the distribution
     *  @return value in distribution that is less than or equal to probability (x)
     */
    double icdf(double x)
    {
        if (x <= 0 || s == 0) {
            return 0;
        }
        return s * log(x / (1.0 - x));
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
