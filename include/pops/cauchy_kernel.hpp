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

#ifndef POPS_CAUCHY_KERNEL_HPP
#define POPS_CAUCHY_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"

#include <random>

namespace pops {

using std::pow;
using std::tan;

/*! Dispersal kernel for log secant
 */
class CauchyKernel
{
protected:
    double s;
    std::cauchy_distribution<double> cauchy_distribution;

public:
    CauchyKernel(double ss, double unused) : s(ss), cauchy_distribution(0, s) {}

    template<class Generator>
    double random(Generator& generator)
    {
        return std::abs(cauchy_distribution(generator));
    }

    double pdf(double x)
    {
        return 1 / ((s * M_PI) * (1 + (pow(x / s, 2))));
    }
    // Inverse cdf (quantile function)
    double icdf(double x)
    {
        return s * tan(M_PI * (x - 0.5));
    }
    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        return type == DispersalKernelType::Cauchy;
    }
};

}  // namespace pops

#endif  // POPS_CAUCHY_KERNEL_HPP
