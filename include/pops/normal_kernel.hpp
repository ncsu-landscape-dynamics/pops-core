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

#ifndef POPS_NORMAL_KERNEL_HPP
#define POPS_NORMAL_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"

#include <random>

namespace pops {

using std::pow;

//  0.147 used for a relative error of about 2*10^-3
//  equation for approximation for erfinv from
//  "A handy approximation for the error function and its inverse" by Sergei Winitzki
float inv_erf(float x)
{
    float sign = (x < 0) ? -1.0f : 1.0f;

    float b = 2 / (M_PI * 0.147) + 0.5f * log(1 - pow(x, 2));

    return (sign * sqrt(-b + sqrt(pow(b, 2) - (1 / (0.147) * log(1 - pow(x, 2))))));
}

/*! Dispersal kernel for log secant
 */
class NormalKernel
{
protected:
    double mu;
    double sigma;
    std::normal_distribution<double> normal_distribution;

public:
    NormalKernel(double m, double s) : mu(m), sigma(s), normal_distribution(sigma, mu)
    {}

    template<class Generator>
    double random(Generator& generator)
    {
        return std::abs(normal_distribution(generator));
    }

    double pdf(double x)
    {
        if (mu == 0 || sigma == 1) {
            return 1 / (sqrt(2 * M_PI)) * exp(-0.5 * pow(x, 2));
        }
        return 1 / (sigma * sqrt(2 * M_PI)) * exp(-0.5 * pow((x - mu) / sigma, 2));
    }

    double icdf(double x)
    {
        if (x <= 0) {
            return 0;
        }
        return mu + (sigma * std::sqrt(2) * (inv_erf((2 * x) - 1)));
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        return type == DispersalKernelType::Normal;
    }
};

}  // namespace pops

#endif  // POPS_NORMAL_KERNEL_HPP