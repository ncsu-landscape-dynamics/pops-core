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

#ifndef POPS_EXPONENTIAL_POWER_KERNEL_HPP
#define POPS_EXPONENTIAL_POWER_KERNEL_HPP

#include "kernel_types.hpp"
#include "deterministic_kernel.hpp"
#include "gamma_kernel.hpp"
#include "raster.hpp"

#include <random>

namespace pops {

using std::pow;

/*! Dispersal kernel for exponential power distribution
 */
class ExponentialPowerKernel
{
protected:
    double alpha;
    double beta;

public:
    ExponentialPowerKernel(double a, double b) : alpha(a), beta(b) {}

    template<class Generator>
    double random(Generator& generator)
    {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        double x = distribution(generator);
        // GammaDistribution gamma_distribution(1.0 / beta_, 1.0 / pow(alpha_, beta_));
        // double gamma = gamma_distribution.icdf(2 * std::abs(x - 0.5));
        // return (x - 0.5) * pow(gamma, 1.0 / beta_);
        return icdf(x);
    }

    double pdf(double x)
    {
        if (beta == 0) {
            return 0;
        }
        return (beta / (2 * alpha * std::tgamma(1.0 / beta)))
               * pow(exp(-x / alpha), beta);
    }

    double icdf(double x)
    {
        GammaDistribution gamma_distribution(1.0 / beta, 1.0 / pow(alpha, beta));
        double gamma = gamma_distribution.icdf(2 * std::abs(x - 0.5));
        return (x - 0.5) * pow(gamma, 1.0 / beta);
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        return type == DispersalKernelType::ExponentialPower;
    }
};

}  // namespace pops

#endif  // POPS_EXPONENTIAL_POWER_KERNEL_HPP