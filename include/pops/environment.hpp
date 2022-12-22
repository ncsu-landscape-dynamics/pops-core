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
 * @brief Type of weather
 *
 * This includes all the ways of how weather coefficient for one step can be obtained.
 */
enum class WeatherType
{
    Deterministic,  ///< Weather is taken from a time series
    Probabilistic,  ///< Weather is generated from a distribution
};

/*! Get a corresponding enum value for a string which represents a weather type.
 *
 * Throws an std::invalid_argument exception if the values was not
 * found or is not supported (which is the same thing).
 */
inline WeatherType weather_type_from_string(const std::string& text)
{
    if (text == "Deterministic" || text == "Deterministic")
        return WeatherType::Deterministic;
    if (text == "probabilistic" || text == "Probabilistic")
        return WeatherType::Probabilistic;
    throw std::invalid_argument(
        "weather_type_from_string: Invalid value '" + text + "' provided");
}

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

    /**
     * @brief Get weather coefficient raster
     *
     * @return Reference to the current weather coefficient raster
     *
     * @throw std::logic_error when coefficient is not set
     */
    const FloatRaster& weather_coefficient() const
    {
        if (!current_weather_coefficient) {
            throw std::logic_error("Weather coefficient used, but not provided");
        }
        return *current_weather_coefficient;
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
