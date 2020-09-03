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

#ifndef POPS_EXPONENTIAL_KERNEL_HPP
#define POPS_EXPONENTIAL_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"

#include <random>

namespace pops {

using std::exp;
using std::log;

/*! Dispersal kernel for Exponential distribution
 */
class ExponentialKernel
{
protected:
    double beta;
    std::exponential_distribution<double> exponential_distribution;

public:
    ExponentialKernel(double b, double unused)
        : beta(b), exponential_distribution(1.0 / beta)
    {}

    template<class Generator>
    double random(Generator& generator)
    {
        return std::abs(exponential_distribution(generator));
    }

    // assumes mu is 0 which is traditionally accepted
    double pdf(double x)
    {
        return (1 / beta) * (exp(-x / beta));
    }
    // Inverse cdf (quantile function)
    double icdf(double x)
    {
        if (beta == 1) {
            return -log(1 - x);
        }
        else {
            return -beta * log(1 - x);
        }
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        return type == DispersalKernelType::Exponential;
    }
};

}  // namespace pops

#endif  // POPS_EXPONENTIAL_KERNEL_HPP