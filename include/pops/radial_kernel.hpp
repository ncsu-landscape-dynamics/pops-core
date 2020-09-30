/*
 * PoPS model - radial disperal kernel
 *
 * Copyright (C) 2015-2020 by the authors.
 *
 * Authors: Vaclav Petras (wenzeslaus gmail com)
 *          Chris Jones (cjones1688 gmail com)
 *          Anna Petrasova (kratochanna gmail com)
 *          Zexi Chen (zchen22 ncsu edu)
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_RADIAL_KERNEL_HPP
#define POPS_RADIAL_KERNEL_HPP

#include "kernel_types.hpp"
#include "hyperbolic_secant_kernel.hpp"
#include "logistic_kernel.hpp"
#include "exponential_power_kernel.hpp"
#include "exponential_kernel.hpp"
#include "cauchy_kernel.hpp"
#include "gamma_kernel.hpp"
#include "lognormal_kernel.hpp"
#include "normal_kernel.hpp"
#include "weibull_kernel.hpp"
#include "power_law_kernel.hpp"

#include <cmath>
#include <map>
#include <tuple>
#include <array>
#include <random>
#include <algorithm>
#include <stdexcept>

// PI is used in the code and M_PI is not guaranteed
// fix it, but prefer the system definition
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef PI
#define PI M_PI
#endif

namespace pops {

// we need to bring sqrt to the namespace
// otherwise only candidates are raster-related
using std::sqrt;

/*! Von Mises Distribution (Circular data distribution)

    mu is the mean angle, expressed in radians between 0 and 2*pi,
    and kappa is the concentration parameter, which must be greater
    than or equal to zero. If kappa is equal to zero, this distribution
    reduces to a uniform random angle over the range 0 to 2*pi.
*/
class von_mises_distribution
{
public:
    von_mises_distribution(double mu, double kappa)
        : mu(mu), kappa(kappa), distribution(0.0, 1.0)
    {}
    template<class Generator>
    double operator()(Generator& generator)
    {
        double a, b, c, f, r, theta, u1, u2, u3, z;

        if (kappa <= 1.e-06)
            return 2 * PI * distribution(generator);

        a = 1.0 + sqrt(1.0 + 4.0 * kappa * kappa);
        b = (a - sqrt(2.0 * a)) / (2.0 * kappa);
        r = (1.0 + b * b) / (2.0 * b);

        while (true) {
            u1 = distribution(generator);
            z = cos(PI * u1);
            f = (1.0 + r * z) / (r + z);
            c = kappa * (r - f);
            u2 = distribution(generator);
            if (u2 <= c * (2.0 - c) || u2 < c * exp(1.0 - c))
                break;
        }

        u3 = distribution(generator);
        if (u3 > 0.5) {
            theta = fmod(mu + acos(f), 2 * PI);
        }
        else {
            theta = fmod(mu - acos(f), 2 * PI);
        }
        return theta;
    }

private:
    double mu;
    double kappa;
    std::uniform_real_distribution<double> distribution;
};

/*! Spread direction
 *
 * Spread, typically wind, direction.
 * Values are in degrees and are used in computations.
 * `None` means that there is no wind.
 */
enum class Direction
{
    N = 0,  //!< North
    NE = 45,  //!< Northeast
    E = 90,  //!< NEast
    SE = 135,  //!< Southeast
    S = 180,  //!< South
    SW = 225,  //!< Southwest
    W = 270,  //!< West
    NW = 315,  //!< Northwest
    None  //!< No direction (non-directional)
};

/*! Get a corresponding enum value for a string which direction.
 *
 * Throws an std::invalid_argument exception if the values was not
 * found or is not supported (which is the same thing).
 */
inline Direction direction_from_string(const std::string& text)
{
    std::map<std::string, Direction> mapping{
        {"N", Direction::N},
        {"NE", Direction::NE},
        {"E", Direction::E},
        {"SE", Direction::SE},
        {"S", Direction::S},
        {"SW", Direction::SW},
        {"W", Direction::W},
        {"NW", Direction::NW},
        {"NONE", Direction::None},
        {"None", Direction::None},
        {"none", Direction::None},
        {"", Direction::None}};
    try {
        return mapping.at(text);
    }
    catch (const std::out_of_range&) {
        throw std::invalid_argument(
            "direction_from_string: Invalid"
            " value '"
            + text + "' provided");
    }
}

/*! Overload which allows to pass C-style string which is nullptr (NULL)
 */
inline Direction direction_from_string(const char* text)
{
    // call the string version
    return direction_from_string(text ? std::string(text) : std::string());
}

/*! Dispersal kernel providing all the radial kernels.
 *
 * We understand a radial kernel to be a kernel which has parameters
 * which translate into a distance and direction.
 *
 * To add new kernel which fits with the other kernels supported by this
 * class, add new member, its initialization from parameters,
 * its implementation in the function call operator, and extend the
 * supports_kernel() function.
 */
template<typename IntegerRaster>
class RadialDispersalKernel
{
protected:
    // the west-east resolution of the pixel
    double east_west_resolution;
    // the north-south resolution of the pixel
    double north_south_resolution;
    DispersalKernelType dispersal_kernel_type_;
    CauchyKernel cauchy_distribution;
    ExponentialKernel exponential_distribution;
    WeibullKernel weibull_distribution;
    NormalKernel normal_distribution;
    LogNormalKernel lognormal_distribution;
    PowerLawKernel power_law_distribution;
    HyperbolicSecantKernel hyperbolic_secant_distribution;
    GammaKernel gamma_distribution;
    ExponentialPowerKernel exponential_power_distribution;
    LogisticKernel logistic_distribution;
    von_mises_distribution von_mises;

public:
    RadialDispersalKernel(
        double ew_res,
        double ns_res,
        DispersalKernelType dispersal_kernel,
        double distance_scale,
        Direction dispersal_direction = Direction::None,
        double dispersal_direction_kappa = 0,
        const IntegerRaster& dispersers = {{0}},
        double dispersal_percentage = 0.99,
        double shape = 1)
        : east_west_resolution(ew_res),
          north_south_resolution(ns_res),
          dispersal_kernel_type_(dispersal_kernel),
          // Here we initialize all distributions,
          // although we won't use all of them.
          cauchy_distribution(distance_scale),
          // When lambda is higher, exponential gives less higher values,
          // so we do multiplicative inverse to behave like cauchy.
          exponential_distribution(distance_scale),
          weibull_distribution(distance_scale, shape),
          normal_distribution(distance_scale),
          lognormal_distribution(distance_scale),
          power_law_distribution(distance_scale, shape),
          hyperbolic_secant_distribution(distance_scale),
          gamma_distribution(distance_scale, shape),
          exponential_power_distribution(distance_scale, shape),
          logistic_distribution(distance_scale),
          // if no wind, then kappa is 0
          // TODO: change these two computations to standalone inline
          // functions (dir to rad and adjust kappa)
          von_mises(
              static_cast<int>(dispersal_direction) * PI / 180,
              dispersal_direction == Direction::None ? 0 : dispersal_direction_kappa)
    {}

    /*! Generates a new position for the spread.
     *
     * The randomness is based on the *generator*. The result may depend
     * on previous calls of this operator (see e.g.
     * `std::cauchy_distribution<RealType>::reset()`).
     * Parameters *row* and *col* are row and column position of the
     * current disperser. The generated position will be relative to it.
     */
    template<typename Generator>
    std::tuple<int, int> operator()(Generator& generator, int row, int col)
    {
        double distance = 0;
        double theta = 0;
        // switch between the supported kernels
        // generate the distance from cauchy distribution or cauchy mixture distribution
        if (dispersal_kernel_type_ == DispersalKernelType::Cauchy) {
            distance = std::abs(cauchy_distribution.random(generator));
        }
        else if (dispersal_kernel_type_ == DispersalKernelType::Exponential) {
            distance = std::abs(exponential_distribution.random(generator));
        }
        else if (dispersal_kernel_type_ == DispersalKernelType::Weibull) {
            distance = std::abs(weibull_distribution.random(generator));
        }
        else if (dispersal_kernel_type_ == DispersalKernelType::Normal) {
            distance = std::abs(normal_distribution.random(generator));
        }
        else if (dispersal_kernel_type_ == DispersalKernelType::LogNormal) {
            distance = std::abs(lognormal_distribution.random(generator));
        }
        else if (dispersal_kernel_type_ == DispersalKernelType::PowerLaw) {
            distance = std::abs(power_law_distribution.random(generator));
        }
        else if (dispersal_kernel_type_ == DispersalKernelType::HyperbolicSecant) {
            distance = std::abs(hyperbolic_secant_distribution.random(generator));
        }
        else if (dispersal_kernel_type_ == DispersalKernelType::Gamma) {
            distance = std::abs(gamma_distribution.random(generator));
        }
        else if (dispersal_kernel_type_ == DispersalKernelType::Logistic) {
            distance = std::abs(exponential_power_distribution.random(generator));
        }
        else if (dispersal_kernel_type_ == DispersalKernelType::Logistic) {
            distance = std::abs(logistic_distribution.random(generator));
        }
        else {
            // TODO: move this to constructor (fail fast)
            // not all allowed kernels will/are supported by this class
            throw std::invalid_argument(
                "RadialDispersalKernel: Unsupported dispersal kernel type");
        }
        theta = von_mises(generator);

        row -= round(distance * cos(theta) / north_south_resolution);
        col += round(distance * sin(theta) / east_west_resolution);

        return std::make_tuple(row, col);
    }

    /*! Returns true if the kernel class support a given kernel type
     *
     * \warning This function is experimental and may be removed or
     * changed at any time.
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        static const std::array<DispersalKernelType, 10> supports = {
            DispersalKernelType::Cauchy,
            DispersalKernelType::Exponential,
            DispersalKernelType::Weibull,
            DispersalKernelType::Normal,
            DispersalKernelType::LogNormal,
            DispersalKernelType::PowerLaw,
            DispersalKernelType::HyperbolicSecant,
            DispersalKernelType::Gamma,
            DispersalKernelType::ExponentialPower,
            DispersalKernelType::Logistic};
        auto it = std::find(supports.cbegin(), supports.cend(), type);
        return it != supports.cend();
    }
};

}  // namespace pops

#endif  // POPS_RADIAL_KERNEL_HPP
