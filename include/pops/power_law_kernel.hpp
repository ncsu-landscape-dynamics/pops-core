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

#ifndef POPS_POWER_LAW_KERNEL_HPP
#define POPS_POWER_LAW_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"

#include <random>

namespace pops {

using std::pow;

/*! Dispersal kernel for power law distribution
 */
class PowerLawKernel
{
protected:
    double xmin;
    double alpha;

public:
    PowerLawKernel(double xm, double a) : xmin(xm), alpha(a) {}

    template<class Generator>
    double random(Generator& generator)
    {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        double x = distribution(generator);
        // return pow(x, (1.0 / (-alpha_ + 1.0))) * xmin_;
        // since distribution is 0 to 1 xmin needs to be greater than 0
        if (xmin <= 0) {
            xmin = 0.01;
        }
        return icdf(x);
    }

    // Should only work with alpha > 1
    // Since power law begins at 1 the distribution the center
    // square in the prob matrix is always 0 - should maybe shift the results?
    double pdf(double x)
    {
        if (x <= 0 || xmin == 0 || alpha <= 1.0) {
            return 0;
        }
        return ((alpha - 1.0) / xmin) * pow(x / xmin, -alpha);
    }

    double icdf(double x)
    {
        if (x <= 0 || xmin == 0 || alpha <= 1.0) {
            return 0;
        }
        return pow(x, (1.0 / (-alpha + 1.0))) * xmin;
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        return type == DispersalKernelType::PowerLaw;
    }
};

}  // namespace pops

#endif  // POPS_POWER_LAW_KERNEL_HPP