/*
 * PoPS model - environment for hosts and pests
 *
 * Copyright (C) 2022 by the authors.
 *
 * Authors: Vaclav Petras (wenzeslaus gmail com)
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_ENVIRONMENT_HPP
#define POPS_ENVIRONMENT_HPP

#include <cmath>
#include <memory>
#include <tuple>
#include <vector>
#include <random>
#include <string>
#include <stdexcept>

#include "utils.hpp"

namespace pops {

/**
 * Encapsulates surrounding environment
 *
 * Currently, only handles weather coefficient for soils. Holds only the current state.
 */
template<typename IntegerRaster, typename FloatRaster, typename RasterIndex = int>
class Environment
{
public:
    Environment() {}

    /**
     * @brief Update the current weather coefficient
     *
     * @param raster Raster with the weather coefficient.
     */
    void update_weather_coefficient(const FloatRaster& raster)
    {
        current_weather_coefficient = &raster;
    }

    template<typename Generator>
    void update_weather_from_distribution(
        const FloatRaster& mean, const FloatRaster& stddev, Generator& generator)
    {
        // probably just pseudo-code, see what works with Rcpp
        stored_weather_coefficient = FloatRaster(mean.rows(), mean.cols());
        // possibly use suitable cells here
        for (RasterIndex i = 0; i < mean.rows(); ++i) {
            for (RasterIndex j = 0; i < mean.rows(); ++i) {
                std::normal_distribution<double> distribution{mean(i, j), stddev(i, j)};
                stored_weather_coefficient(i, j) = distribution(generator);
            }
        }
        current_weather_coefficient = &stored_weather_coefficient;
    }

    /**
     * @brief Get weather coefficient at a given cell
     *
     * @param row Cell row number
     * @param col Cell column number
     * @return Current value at the given cell
     *
     * @throw std::logic_error when coefficient is not set
     */
    double weather_coefficient_at(RasterIndex row, RasterIndex col) const
    {
        if (!current_weather_coefficient) {
            throw std::logic_error("Weather coefficient used, but not provided");
        }
        return current_weather_coefficient->operator()(row, col);
    }

protected:
    /**
     * Current weather coefficient
     *
     * Value may not be set and these cases should produce exceptions.
     */
    const FloatRaster* current_weather_coefficient{nullptr};
    FloatRaster stored_weather_coefficient;
};

}  // namespace pops

#endif  // POPS_ENVIRONMENT_HPP
