/*
 * PoPS model - deterministic dispersal kernel
 *
 * Copyright (C) 2015-2020 by the authors.
 *
 * Authors: Margaret Lawrimore (malawrim ncsu edu)
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_DETERMINISTIC_KERNEL_HPP
#define POPS_DETERMINISTIC_KERNEL_HPP

#include <vector>
#include <tuple>

#include "raster.hpp"
#include "kernel_types.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef PI
#define PI M_PI
#endif

namespace pops {

using std::pow;
using std::tan;
using std::exp;
using std::log;
using std::ceil;
using std::abs;
using std::sqrt;

//  0.147 used for a relative error of about 2*10^-3
//  equation for approximation for erfinv from
//  "A handy approximation for the error function and its inverse" by Sergei Winitzki
float inv_erf(float x)
{
    float sign = (x < 0) ? -1.0f : 1.0f;

    float b = 2 / (M_PI * 0.147) + 0.5f * log(1 - pow(x, 2));

    return (sign * sqrt(-b + sqrt(pow(b, 2) - (1 / (0.147) * log(1 - pow(x, 2))))));
}

/*!
 * Cauchy distribution
 * Includes probability density function and inverse cumulative distribution function
 * pdf returns the probability that the variate has the value x
 * icdf returns the upper range that encompasses x percent of the distribution (e.g for
 * 99% input .99)
 */
class CauchyDistribution
{
public:
    CauchyDistribution(double scale) : s(scale) {}

    double pdf(double x)
    {
        return 1 / ((s * M_PI) * (1 + (pow(x / s, 2))));
    }
    // Inverse cdf (quantile function)
    double icdf(double x)
    {
        return s * tan(M_PI * (x - 0.5));
    }

private:
    // scale parameter - 1 for standard
    double s;
};

/*!
 * Exponential distribution
 * Includes probability density function and inverse cumulative distribution function
 * pdf returns the probability that the variate has the value x
 * icdf returns the upper range that encompasses x percent of the distribution (e.g for
 * 99% input 0.99)
 */
class ExponentialDistribution
{
public:
    ExponentialDistribution(double scale) : beta(scale) {}
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

private:
    // scale parameter - 1 for standard
    // equal to 1/lambda
    double beta;
};
/*!
 * Log-Normal distribution
 * Includes probability density function and cumulative distribution function
 * pdf returns the probability that the variate has the value x
 * cdf returns the upper range that encompasses x percent of the distribution (e.g
 * for 99% input 0.99)
 */
class LogNormalDistribution
{
public:
    LogNormalDistribution(double scale, double locator) : mu(scale), sigma(locator) {}

    double pdf(double x)
    {
        if (x <= 0 || sigma == 0) {
            return 0;
        }
        return (1 / x) * (1 / (sigma * sqrt(2 * M_PI)))
               * exp(-(pow(log(x) - mu, 2)) / (2 * pow(sigma, 2)));
    }

    double icdf(double x)
    {
        if (x <= 0) {
            return 0;
        }
        return exp(mu + sqrt(2 * pow(sigma, 2)) * inv_erf((2 * x) - 1));
        // double z = 0.0;
        // double count = 0.0;
        // while ( z < x ) {
        //    count += 0.5;
        //    z = 0.5 * std::erfc(-(log(count) - mu) / (sigma * sqrt(2)));
        //}
        // return count;
    }

private:
    double sigma;  // standard deviation
    double mu;  // mean
};
/*!
 * Gamma distribution
 * Includes probability density function and cumulative distribution function
 * pdf returns the probability that the variate has the value x
 * cdf returns the upper range that encompasses x percent of the distribution (e.g for
 * 99% input 0.99)
 */
class GammaDistribution
{
public:
    GammaDistribution(double scale, double locator) : alpha(scale), theta(locator) {}

    double pdf(double x)
    {
        if (x < 0 || alpha < 0 || theta < 0) {
            return 0;
        }
        return 1.0 / (std::tgamma(alpha) * pow(theta, alpha)) * pow(x, (alpha - 1))
               * exp(-x / theta);
    }

    double cdf(double x)
    {
        double sum = 0.0;
        double beta = 1.0 / theta;
        for (int i = 0; i < alpha; i++) {
            // tgamma = (i-1)! used since c++ has no factorial in std lib
            sum += pow(beta * x, i) / std::tgamma(i + 1);
        }
        return 1 - sum * exp(-beta * x);
    }

    // there is no known closed-form solution
    // R and MatLab use an iterative approach (Newton's method)
    double icdf(double x)
    {
        // pick starting approximation using lognormal icdf
        LogNormalDistribution lognormal(0, 1);
        double guess = lognormal.icdf(x);
        // TODO add cdf function
        double check = cdf(guess);
        double numiterations = 500;  // will need to adjust this
        double precision = 0.001;  // will need to adjust this
        for (int i = 0; i < numiterations; i++) {
            if (check < (x - precision) || check > (x + precision)) {
                double dif = check - x;
                // if dif is positive guess is greater than needed
                // if dif is negative guess is less than needed
                double past_guess = guess;
                double derivative = (check - x) / pdf(guess);
                // limit size of next guess
                guess = std::max(guess / 10, std::min(guess * 10, guess - derivative));
                check = cdf(guess);
                // Check if we went to far and need to backtrack
                for (int j = 0; j < 10; j++) {
                    if (std::abs(dif) < std::abs(check - x)) {
                        past_guess = guess;
                        guess = (guess + past_guess) / 2.0;
                        check = cdf(guess);
                    }
                }
            }
            else {
                return guess;
            }
        }
        // TODO error message
        std::cout << "unable to find solution to gamma icdf(" << x << ")\n";
        return -1;
    }

private:
    double alpha;  // scale
    double theta;  // shape
};

/*!
 * Exponential Power distribution
 * Includes probability density function and cumulative distribution function
 * pdf returns the probability that the variate has the value x
 * icdf returns the upper range that encompasses x percent of the distribution (e.g
 * for 99% input 0.99)
 */
class ExponentialPowerDistribution
{
public:
    ExponentialPowerDistribution(double scale, double locator)
        : alpha(scale), beta(locator)
    {}

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

private:
    double alpha;  // scale
    double beta;  // shape
};
/*!
 * Weibull distribution
 * Includes probability density function and cumulative distribution function
 * pdf returns the probability that the variate has the value x
 * icdf returns the upper range that encompasses x percent of the distribution (e.g
 * for 99% input 0.99)
 */
class WeibullDistribution
{
public:
    WeibullDistribution(double scale, double locator) : b(scale), a(locator) {}

    double pdf(double x)
    {
        if (x < 0 || b <= 0 || a < 0) {
            return 0;
        }
        // Note: If the value inside exp() is too large it returns zero
        // any a >= 2 returns zero
        return ((a / b) * pow(x / b, a - 1) * exp(-pow(x / b, a)));
    }

    double icdf(double x)
    {
        if (x < 0 || x >= 1 || b <= 0 || a < 0) {
            return 0;
        }
        return a * pow(-(log(1 - x)), (1.0 / b));
    }

private:
    double a;  // shape
    double b;  // scale
};

/*!
 * Normal distribution
 * Includes probability density function and cumulative distribution function
 * pdf returns the probability that the variate has the value x
 * icdf returns the upper range that encompasses x percent of the distribution (e.g
 * for 99% input 0.99)
 */
class NormalDistribution
{
public:
    NormalDistribution(double scale, double locator) : mu(scale), sigma(locator) {}

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
        // double z = 0.0;
        // double count = 0.0;
        // while ( z < x ) {
        //    count += 0.5;
        //    z = 0.5 * (1 + std::erf((count - mu) / (sigma * sqrt(2))));
        //}
        // return count;
    }

private:
    double sigma;  // standard deviation
    double mu;  // mean
};

/*!
 * Power law distribution
 * Includes probability density function and cumulative distribution function
 * pdf returns the probability that the variate has the value x
 * cdf returns the upper range that encompasses x percent of the distribution (e.g
 * for 99% input 0.99)
 */
class PowerLawDistribution
{
public:
    PowerLawDistribution(double scale, double locator) : alpha(scale), xmin(locator) {}

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

private:
    double alpha;
    double xmin;
};

// Note: Can't find anything called log-hyperbolic secant
/*!
 * Hyperbolic Secant distribution
 * Includes probability density function and cumulative distribution function
 * pdf returns the probability that the variate has the value x
 * cdf returns the upper range that encompasses x percent of the distribution (e.g
 * for 99% input 0.99)
 */
class HyperbolicSecantDistribution
{
public:
    HyperbolicSecantDistribution(double scale, double locator)
        : mu(scale), sigma(locator)
    {}

    double pdf(double x)
    {
        if (x <= 0 || sigma == 0) {
            return 0;
        }
        if (mu == 0 && sigma == 1) {
            return 0.5 * (1 / cosh((M_PI * x) / 2));
        }
        return (1 / (2 * sigma)) * (1 / cosh((M_PI * (x - mu)) / (2 * sigma)));
    }

    double icdf(double x)
    {
        if (x <= 0 || sigma == 0) {
            return 0;
        }
        if (mu == 0 && sigma == 1) {
            return (2 / M_PI) * log(tan(M_PI / 2 * x));
        }
        return ((log(tan((x * M_PI) / 2)) * (2 * sigma)) / M_PI) + mu;
    }

private:
    double mu;  //
    double sigma;  //
};

/*!
 * Logistic distribution
 * Includes probability density function and cumulative distribution function
 * pdf returns the probability that the variate has the value x
 * cdf returns the upper range that encompasses x percent of the distribution (e.g
 * for 99% input 0.99)
 */
class LogisticDistribution
{
public:
    LogisticDistribution(double scale, double locator) : mu(scale), s(locator) {}

    double pdf(double x)
    {
        if (x < 0 || s == 0) {
            return 0;
        }
        if (mu == 0 && s == 1) {
            return exp(-x) / pow(1 + exp(-x), 2);
        }
        return (exp(-(x - mu) / s)) / (s * pow(1 + exp(-(x - mu) / s), 2));
    }

    double icdf(double x)
    {
        if (x <= 0 || s == 0) {
            return 0;
        }
        if (mu == 0 && s == 1) {
            return log(x / (1 - x));
        }
        return mu + (s * log(x / (1 - x)));
    }

private:
    double mu;  //
    double s;  //
};

/*!
 * Dispersal kernel for deterministic spread to cell with highest probability of
 * spread
 *
 * Dispersal Kernel type determines use of Exponential or Cauchy distribution
 * to find probability.
 *
 * dispersal_percentage is the percent of all possible dispersal to be included
 * in the moving window size (e.g for 99% input 0.99).
 *
 * Useful for testing as it is deterministic and provides fully replicable results
 */
template<typename IntegerRaster>
class DeterministicDispersalKernel
{
protected:
    const IntegerRaster& dispersers_;
    // row/col position of middle cell
    int mid_row = 0;
    int mid_col = 0;
    // position of cell from previous call
    int prev_row = -1;
    int prev_col = -1;
    // number of rows/cols in the probability window
    int number_of_rows = 0;
    int number_of_columns = 0;
    // maximum distance from center cell to outer cells
    double max_distance{0};
    Raster<double> probability;
    Raster<double> probability_copy;
    CauchyDistribution cauchy;
    ExponentialDistribution exponential;
    WeibullDistribution weibull;
    LogNormalDistribution logNormal;
    NormalDistribution normal;
    HyperbolicSecantDistribution hyperbolicSecant;
    PowerLawDistribution powerLaw;
    GammaDistribution gamma;
    ExponentialPowerDistribution exponentialPower;
    LogisticDistribution logistic;

    DispersalKernelType kernel_type_;
    double proportion_of_dispersers;
    // the west-east resolution of the pixel
    double east_west_resolution;
    // the north-south resolution of the pixel
    double north_south_resolution;

public:
    DeterministicDispersalKernel(
        DispersalKernelType dispersal_kernel,
        const IntegerRaster& dispersers,
        double dispersal_percentage,
        double ew_res,
        double ns_res,
        double distance_scale,
        double locator)
        : dispersers_(dispersers),
          cauchy(distance_scale),
          exponential(distance_scale),
          weibull(distance_scale, locator),
          logNormal(distance_scale, locator),
          normal(distance_scale, locator),
          hyperbolicSecant(distance_scale, locator),
          powerLaw(distance_scale, locator),
          logistic(distance_scale, locator),
          gamma(distance_scale, locator),
          exponentialPower(distance_scale, locator),
          kernel_type_(dispersal_kernel),
          east_west_resolution(ew_res),
          north_south_resolution(ns_res)
    {
        // We initialize max distance only for the supported kernels.
        // For the others, we report the error only when really called
        // to allow use of this class in initialization phase.
        if (kernel_type_ == DispersalKernelType::Cauchy) {
            max_distance = cauchy.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::Exponential) {
            max_distance = exponential.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::Weibull) {
            max_distance = weibull.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::Normal) {
            max_distance = normal.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::LogNormal) {
            max_distance = logNormal.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::HyperbolicSecant) {
            max_distance = hyperbolicSecant.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::PowerLaw) {
            max_distance = powerLaw.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::Logistic) {
            max_distance = logistic.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::Gamma) {
            max_distance = gamma.icdf(dispersal_percentage);
        }
        else if (kernel_type_ == DispersalKernelType::ExponentialPower) {
            max_distance = exponentialPower.icdf(dispersal_percentage);
        }
        number_of_columns = ceil(max_distance / east_west_resolution) * 2 + 1;
        number_of_rows = ceil(max_distance / north_south_resolution) * 2 + 1;
        Raster<double> prob_size(number_of_rows, number_of_columns, 0);
        probability = prob_size;
        probability_copy = prob_size;
        mid_row = number_of_rows / 2;
        mid_col = number_of_columns / 2;
        double sum = 0.0;
        for (int i = 0; i < number_of_rows; i++) {
            for (int j = 0; j < number_of_columns; j++) {
                double distance_to_center = std::sqrt(
                    pow((abs(mid_row - i) * east_west_resolution), 2)
                    + pow((abs(mid_col - j) * north_south_resolution), 2));
                // determine probability based on distance
                if (kernel_type_ == DispersalKernelType::Cauchy) {
                    probability(i, j) = abs(cauchy.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::Exponential) {
                    probability(i, j) = abs(exponential.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::Weibull) {
                    probability(i, j) = abs(weibull.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::Normal) {
                    probability(i, j) = abs(normal.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::LogNormal) {
                    probability(i, j) = abs(logNormal.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::PowerLaw) {
                    probability(i, j) = abs(powerLaw.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::HyperbolicSecant) {
                    probability(i, j) = abs(hyperbolicSecant.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::Logistic) {
                    probability(i, j) = abs(logistic.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::Gamma) {
                    probability(i, j) = abs(gamma.pdf(distance_to_center));
                }
                else if (kernel_type_ == DispersalKernelType::ExponentialPower) {
                    probability(i, j) = abs(exponentialPower.pdf(distance_to_center));
                }
                sum += probability(i, j);
            }
        }
        // normalize based on the sum of all probabilities in the raster
        probability /= sum;
    }

    /*! Generates a new position for the spread.
     *
     *  Creates a copy of the probability matrix to mark where dispersers are
     * assigned. New window created any time a new cell is selected from
     * simulation.disperse
     *
     *  Selects next row/col value based on the cell with the highest probability
     *  in the window.
     *
     */
    template<class Generator>
    std::tuple<int, int> operator()(Generator& generator, int row, int col)
    {
        if (kernel_type_ != DispersalKernelType::Cauchy
            && kernel_type_ != DispersalKernelType::Exponential
            && kernel_type_ != DispersalKernelType::Weibull
            && kernel_type_ != DispersalKernelType::Normal
            && kernel_type_ != DispersalKernelType::LogNormal
            && kernel_type_ != DispersalKernelType::HyperbolicSecant
            && kernel_type_ != DispersalKernelType::PowerLaw
            && kernel_type_ != DispersalKernelType::Logistic
            && kernel_type_ != DispersalKernelType::Gamma
            && kernel_type_ != DispersalKernelType::ExponentialPower) {
            throw std::invalid_argument(
                "DeterministicDispersalKernel: Unsupported dispersal kernel type");
        }
        // reset the window if considering a new cell
        if (row != prev_row || col != prev_col) {
            proportion_of_dispersers = 1.0 / (double)dispersers_(row, col);
            probability_copy = probability;
        }

        int row_movement = 0;
        int col_movement = 0;

        double max = (double)-std::numeric_limits<int>::max();
        int max_prob_row = 0;
        int max_prob_col = 0;

        // find cell with highest probability
        for (int i = 0; i < number_of_rows; i++) {
            for (int j = 0; j < number_of_columns; j++) {
                if (probability_copy(i, j) > max) {
                    max = probability_copy(i, j);
                    max_prob_row = i;
                    max_prob_col = j;
                    row_movement = i - mid_row;
                    col_movement = j - mid_col;
                }
            }
        }

        // subtracting 1/number of dispersers ensures we always move the same
        // proportion of the individuals to each cell no matter how many are
        // dispersing
        probability_copy(max_prob_row, max_prob_col) -= proportion_of_dispersers;
        prev_row = row;
        prev_col = col;

        // return values in terms of actual location
        return std::make_tuple(row + row_movement, col + col_movement);
    }
};

}  // namespace pops

#endif  // POPS_DETERMINISTIC_KERNEL_HPP
