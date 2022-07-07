/*
 * PoPS model - soils pool
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

#ifndef POPS_SOILS_HPP
#define POPS_SOILS_HPP

#include <cmath>
#include <memory>
#include <tuple>
#include <vector>
#include <random>
#include <string>
#include <stdexcept>

#include "utils.hpp"
#include "environment.hpp"

namespace pops {

template<typename IntegerRaster, typename FloatRaster, typename RasterIndex = int>
class SoilPool
{
public:
    SoilPool(
        std::vector<IntegerRaster>& rasters,
        const Environment<IntegerRaster, FloatRaster, RasterIndex>& environment,
        bool generate_stochasticity = true,
        bool establishment_stochasticity = true,
        double fixed_soil_probability = 0)
        : rasters_(&rasters),
          environment_(&environment),
          generate_stochasticity_(generate_stochasticity),
          establishment_stochasticity_(establishment_stochasticity),
          fixed_soil_probability_(fixed_soil_probability)
    {
        if (rasters.empty()) {
            throw std::logic_error(
                "List of rasters of SpoilPool needs to have at least one item");
        }
    }

    template<typename Generator>
    int dispersers_from(RasterIndex row, RasterIndex col, Generator& generator)
    {
        auto count = this->total_at(row, col);
        double lambda = environment_->weather_coefficient_at(row, col);
        if (this->generate_stochasticity_) {
            std::poisson_distribution<int> distribution(lambda);
            int dispersers = 0;
            for (int k = 0; k < count; k++) {
                dispersers += distribution(generator);
            }
            return dispersers;
        }
        else {
            return lambda * count;
        }
        // TODO: subtract as with mortality moves
    }

    template<typename Generator>
    void disperser_to(RasterIndex row, RasterIndex col, Generator& generator)
    {
        double current_probability = environment_->weather_coefficient_at(row, col);
        double tester = 1 - fixed_soil_probability_;
        if (this->establishment_stochasticity_)
            tester = distribution_uniform_(generator);
        if (tester < current_probability) {
            this->add_at(row, col);
        }
    }

    template<typename Generator>
    void dispersers_to(
        int dispersers, RasterIndex row, RasterIndex col, Generator& generator)
    {
        for (int i = 0; i < dispersers; i++)
            this->disperser_to(row, col, generator);
    }

    void add_at(RasterIndex row, RasterIndex col, int value = 1)
    {
        rasters_->back()(row, col) += value;
    }

    int total_at(RasterIndex row, RasterIndex col)
    {
        int total = 0;
        for (auto& raster : *rasters_) {
            total += raster(row, col);
        }
        return total;
    }

    void next_step(int step)
    {
        UNUSED(step);
        rotate_left_by_one(*rasters_);
        rasters_->back().fill(0);
    }

protected:
    std::vector<IntegerRaster>* rasters_{nullptr};
    const Environment<IntegerRaster, FloatRaster, RasterIndex>* environment_{nullptr};
    bool generate_stochasticity_{false};
    bool establishment_stochasticity_{false};
    double fixed_soil_probability_{0};
    std::uniform_real_distribution<double> distribution_uniform_{0.0, 1.0};
};

}  // namespace pops

#endif  // POPS_SOILS_HPP
