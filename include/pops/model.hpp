/*
 * Tests for the PoPS Model class.
 *
 * Copyright (C) 2020 by the authors.
 *
 * Authors: Vaclav Petras <wenzeslaus gmail com>
 *
 * This file is part of PoPS.

 * PoPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * PoPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with PoPS. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef POPS_MODEL_HPP
#define POPS_MODEL_HPP

#include "config.hpp"
#include "treatments.hpp"
#include "spread_rate.hpp"
#include "simulation.hpp"
#include "kernel.hpp"
#include "scheduling.hpp"

#include <vector>

namespace pops {

template<typename IntegerRaster, typename FloatRaster, typename RasterIndex>
class Model
{
private:
    Config config_;
    DispersalKernelType natural_kernel;
    DispersalKernelType anthro_kernel;
    UniformDispersalKernel uniform_kernel;
    DeterministicNeighborDispersalKernel natural_neighbor_kernel;
    DeterministicNeighborDispersalKernel anthro_neighbor_kernel;
    Simulation<IntegerRaster, FloatRaster, RasterIndex> simulation_;
    unsigned last_index{0};

public:
    Model(const Config& config)
        : config_(config),
          natural_kernel(kernel_type_from_string(config.natural_kernel_type)),
          anthro_kernel(kernel_type_from_string(config.anthro_kernel_type)),
          uniform_kernel(config.rows, config.cols),
          natural_neighbor_kernel(direction_from_string(config.natural_direction)),
          anthro_neighbor_kernel(direction_from_string(config.anthro_direction)),
          simulation_(
              config.random_seed,
              config.rows,
              config.cols,
              model_type_from_string(config.model_type),
              config.latency_period_steps,
              config.generate_stochasticity,
              config.establishment_stochasticity,
              config.movement_stochasticity)
    {}

    void run_step(
        int step,
        int weather_step,  // TODO: this should be schedule (?)
        IntegerRaster& infected,
        IntegerRaster& susceptible,
        IntegerRaster& total_hosts,  // TODO: How it is with updating this?
        IntegerRaster& dispersers,
        std::vector<IntegerRaster>& exposed,
        std::vector<IntegerRaster>& mortality_tracker,
        IntegerRaster& died,
        const std::vector<FloatRaster>& temperatures,
        const std::vector<FloatRaster>& weather_coefficients,
        Treatments<IntegerRaster, FloatRaster>& treatments,
        IntegerRaster& resistant,
        std::vector<std::tuple<int, int>>& outside_dispersers,  // out
        SpreadRate<IntegerRaster>& spread_rate,  // out
        std::vector<std::vector<int>> movements)
    {
        RadialDispersalKernel<IntegerRaster> natural_radial_kernel(
            config_.ew_res,
            config_.ns_res,
            natural_kernel,
            config_.natural_scale,
            direction_from_string(config_.natural_direction),
            config_.natural_kappa,
            config_.deterministic,
            dispersers,
            config_.dispersal_percentage);
        RadialDispersalKernel<IntegerRaster> anthro_radial_kernel(
            config_.ew_res,
            config_.ns_res,
            anthro_kernel,
            config_.anthro_scale,
            direction_from_string(config_.anthro_direction),
            config_.anthro_kappa,
            config_.deterministic,
            dispersers,
            config_.dispersal_percentage);
        SwitchDispersalKernel<IntegerRaster> natural_selectable_kernel(
            natural_kernel,
            natural_radial_kernel,
            uniform_kernel,
            natural_neighbor_kernel);
        SwitchDispersalKernel<IntegerRaster> anthro_selectable_kernel(
            anthro_kernel,
            anthro_radial_kernel,
            uniform_kernel,
            anthro_neighbor_kernel);
        DispersalKernel<IntegerRaster> dispersal_kernel(
            natural_selectable_kernel,
            anthro_selectable_kernel,
            config_.use_anthropogenic_kernel,
            config_.percent_natural_dispersal);
        int mortality_simulation_year =
            simulation_step_to_action_step(config_.mortality_schedule(), step);
        // removal of dispersers due to lethal tempearture
        if (config_.use_lethal_temperature && config_.lethal_schedule()[step]) {
            int lethal_step =
                simulation_step_to_action_step(config_.lethal_schedule(), step);
            simulation_.remove(
                infected,
                susceptible,
                temperatures[lethal_step],
                config_.lethal_temperature);
        }
        // actual spread
        if (config_.spread_schedule()[step]) {
            simulation_.generate(
                dispersers,
                infected,
                config_.weather,
                weather_coefficients[weather_step],
                config_.reproductive_rate);

            simulation_.disperse_and_infect(
                step,
                dispersers,
                susceptible,
                exposed,
                infected,
                mortality_tracker[mortality_simulation_year],
                total_hosts,
                outside_dispersers,
                config_.weather,
                weather_coefficients[weather_step],
                dispersal_kernel,
                config_.establishment_probability);
            if (config_.use_movements) {
                last_index = simulation_.movement(
                    infected,
                    susceptible,
                    mortality_tracker[mortality_simulation_year],
                    total_hosts,
                    step,
                    last_index,
                    movements,
                    config_.movement_schedule);
            }
        }
        // treatments
        if (config_.use_treatments) {
            bool managed = treatments.manage(step, infected, susceptible, resistant);
            if (managed && config_.use_mortality) {
                // same conditions as the mortality code below
                // TODO: make the mortality timing available as a separate function in
                // the library or simply go over all valid cohorts
                if (mortality_simulation_year >= config_.first_mortality_year - 1) {
                    auto max_index =
                        mortality_simulation_year - (config_.first_mortality_year - 1);
                    for (int age = 0; age <= max_index; age++) {
                        treatments.manage_mortality(step, mortality_tracker[age]);
                    }
                }
            }
        }
        if (config_.use_mortality && config_.mortality_schedule()[step]) {
            // only run to the current year of simulation
            // (first year is 0):
            //   max index == sim year
            // reduced by first time when trees start dying
            // (counted from 1: first year == 1)
            // e.g. for sim year 3, year dying 4, max index is 0
            //   max index = sim year - (dying year - 1)
            // index is negative before we reach the year
            // (so we can skip these years)
            // sim year - (dying year - 1) < 0
            // sim year < dying year - 1
            // TODO: died.zero(); should be done by the caller if needed, document!
            simulation_.mortality(
                infected,
                config_.mortality_rate,
                mortality_simulation_year,
                config_.first_mortality_year - 1,
                died,
                mortality_tracker);
        }
        // compute spread rate
        if (config_.spread_rate_schedule()[step]) {
            unsigned simulation_year =
                simulation_step_to_action_step(config_.spread_rate_schedule(), step);
            spread_rate.compute_yearly_spread_rate(infected, simulation_year);
        }
    }
};

}  // namespace pops

#endif  // POPS_MODEL_HPP
