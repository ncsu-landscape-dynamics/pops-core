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

#include "normal_distribution_with_uniform_fallback.hpp"
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
    None,  ///< No weather
};

/*! Get a corresponding enum value for a string which represents a weather type.
 *
 * Throws an std::invalid_argument exception if the values was not
 * found or is not supported (which is the same thing).
 */
inline WeatherType weather_type_from_string(const std::string& text)
{
    if (text == "deterministic" || text == "Deterministic")
        return WeatherType::Deterministic;
    if (text == "probabilistic" || text == "Probabilistic")
        return WeatherType::Probabilistic;
    if (text.empty() || text == "none" || text == "None" || text == "NONE")
        return WeatherType::None;
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
        if (mean.rows() != stddev.rows()) {
            throw std::invalid_argument(
                "Mean and stddev need to have the same number of rows ("
                + std::to_string(mean.rows()) + " != " + std::to_string(stddev.rows())
                + ")");
        }
        if (mean.cols() != stddev.cols()) {
            throw std::invalid_argument(
                "Mean and stddev need to have the same number of columns ("
                + std::to_string(mean.cols()) + " != " + std::to_string(stddev.cols())
                + ")");
        }
        stored_weather_coefficient = FloatRaster(mean.rows(), mean.cols());
        // possibly use suitable cells here
        for (RasterIndex i = 0; i < mean.rows(); ++i) {
            for (RasterIndex j = 0; j < mean.cols(); ++j) {
                auto mean_value = mean(i, j);
                // In general, to get a specific shape of the distribution, mean can be
                // anything, but we limit that assuming that mean which is out of the
                // desired coefficient range is an erroneous value.
                // Notably, this test is not perfomed for deterministic weather.
                if (mean_value < weather_coefficient_min
                    || mean_value > weather_coefficient_max) {
                    throw std::invalid_argument(
                        std::string("Weather coefficient mean is expected to be ")
                        + "between " + std::to_string(weather_coefficient_min) + " and "
                        + std::to_string(weather_coefficient_max) + ", but is "
                        + std::to_string(mean_value) + " at (" + std::to_string(i)
                        + ", " + std::to_string(j) + ")");
                }
                NormalDistributionWithUniformFallback<double> distribution{
                    mean_value,
                    stddev(i, j),
                    weather_coefficient_min,
                    weather_coefficient_max};
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
    static constexpr double weather_coefficient_min = 0;
    static constexpr double weather_coefficient_max = 1;

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
