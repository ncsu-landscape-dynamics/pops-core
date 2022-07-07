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

template<typename IntegerRaster, typename FloatRaster, typename RasterIndex = int>
class Environment
{
public:
    Environment() {}

    void update_weather_coefficient(const FloatRaster& raster)
    {
        current_weather_coefficient = &raster;
    }

    double weather_coefficient_at(RasterIndex row, RasterIndex col) const
    {
        if (!current_weather_coefficient) {
            throw std::logic_error("Weather coefficient used, but not provided");
        }
        return current_weather_coefficient->operator()(row, col);
    }

protected:
    const FloatRaster* current_weather_coefficient{nullptr};
};

}  // namespace pops

#endif  // POPS_ENVIRONMENT_HPP
