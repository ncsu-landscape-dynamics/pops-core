/*
 * Tests for the PoPS Model class.
 *
 * Copyright (C) 2020-2021 by the authors.
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
#include "switch_kernel.hpp"
#include "kernel.hpp"
#include "scheduling.hpp"
#include "quarantine.hpp"

#include <vector>

namespace pops {

template<
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename Generator = std::default_random_engine,
    typename KernelFactory = DispersalKernel<Generator>(
        const Config&, const IntegerRaster&, const Network<RasterIndex>&)>
class Model
{
protected:
    Config config_;
    DispersalKernelType natural_kernel;
    DispersalKernelType anthro_kernel;
    UniformDispersalKernel uniform_kernel;
    DeterministicNeighborDispersalKernel natural_neighbor_kernel;
    DeterministicNeighborDispersalKernel anthro_neighbor_kernel;
    Simulation<IntegerRaster, FloatRaster, RasterIndex, Generator> simulation_;
    KernelFactory& kernel_factory_;
    unsigned last_index{0};

    /**
     * @brief Create overpopulation movement kernel
     *
     * Same as the natural kernel. The natural kernel parameters are used,
     * but the scale for radial and deterministic kernel is multiplied
     * by the leaving scale coefficient.
     *
     * @param dispersers The disperser raster (reference, for deterministic kernel)
     * @param network Network (initialized or not)
     * @return Created kernel
     */
    SwitchDispersalKernel<IntegerRaster, RasterIndex>
    create_overpopulation_movement_kernel(
        const IntegerRaster& dispersers, const Network<RasterIndex>& network)
    {
        RadialDispersalKernel<IntegerRaster> radial_kernel(
            config_.ew_res,
            config_.ns_res,
            natural_kernel,
            config_.natural_scale * config_.leaving_scale_coefficient,
            direction_from_string(config_.natural_direction),
            config_.natural_kappa,
            config_.shape);
        DeterministicDispersalKernel<IntegerRaster> deterministic_kernel(
            natural_kernel,
            dispersers,
            config_.dispersal_percentage,
            config_.ew_res,
            config_.ns_res,
            config_.natural_scale * config_.leaving_scale_coefficient,
            config_.shape);
        NetworkDispersalKernel<RasterIndex> network_kernel(
            network, config_.network_min_time, config_.network_max_time);
        SwitchDispersalKernel<IntegerRaster, RasterIndex> selectable_kernel(
            natural_kernel,
            radial_kernel,
            deterministic_kernel,
            uniform_kernel,
            network_kernel,
            natural_neighbor_kernel,
            config_.deterministic);
        return selectable_kernel;
    }

public:
    Model(
        const Config& config,
        KernelFactory& kernel_factory =
            create_dynamic_kernel<Generator, IntegerRaster, RasterIndex>)
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
              config.movement_stochasticity),
          kernel_factory_(kernel_factory)
    {}

    /**
     * @brief Run one step of the simulation.
     *
     * The *total_populations* can be total number of hosts in the basic case
     * or it can be the total size of population of all relevant species
     * both host and non-host if dilution effect should be applied.
     * When movement is applied, *total_populations* needs to be only the
     * total number of hosts because movement does not support non-host
     * individuals.
     *
     * No treatment can be applied when movement is active because host movement does
     * not support resistant hosts.
     *
     * *dispersers* is for internal use and for tracking dispersers creation.
     * The input values are ignored and the output is not the current existing
     * dispersers, but only the number of dispersers generated (and subsequently
     * used) in this step. There are no dispersers in between simulation steps.
     *
     * @param step Step number in the simulation.
     * @param[in,out] infected Infected hosts
     * @param[in,out] susceptible Susceptible hosts
     * @param[in,out] total_hosts All host individuals in the area. Is equal to
     * infected + exposed + susceptible in the cell.
     * @param[in,out] total_populations All host and non-host individuals in the area
     * @param[out] dispersers Dispersing individuals (used internally)
     * @param[in,out] total_exposed Sum of all exposed hosts (if SEI model is active)
     * @param[in,out] exposed Exposed hosts (if SEI model is active)
     * @param[in,out] mortality_tracker Mortality tracker used to generate *died*.
     * Expectation is that mortality tracker is of length (1/mortality_rate +
     * mortality_time_lag)
     * @param[out] died Infected hosts which died this step based on the mortality
     * schedule
     * @param[in] temperatures Vector of temperatures used to evaluate lethal
     * temperature
     * @param[in] weather_coefficient Weather coefficient (for the current step)
     * @param[in,out] treatments Treatments to be applied (also tracks use of
     * treatments)
     * @param[in,out] resistant Resistant hosts (host temporarily removed from
     * susceptible hosts)
     * @param[in,out] outside_dispersers Dispersers escaping the rasters (adds to the
     * vector)
     * @param[in,out] spread_rate Spread rate tracker
     * @param[in,out] quarantine Quarantine escape tracker
     * @param[in] quarantine_areas Quarantine areas
     * @param[in] movements Table of host movements
     * @param network Network (initialized or Network::null_network() if unused)
     * @param[in,out] suitable_cells List of indices of cells with hosts
     *
     * @note The parameters roughly correspond to Simulation::disperse()
     * and Simulation::disperse_and_infect() functions, so these can be used
     * for further reference.
     */
    void run_step(
        int step,
        IntegerRaster& infected,
        IntegerRaster& susceptible,
        IntegerRaster& total_populations,
        IntegerRaster& total_hosts,
        IntegerRaster& dispersers,
        IntegerRaster& total_exposed,
        std::vector<IntegerRaster>& exposed,
        std::vector<IntegerRaster>& mortality_tracker,
        IntegerRaster& died,
        const std::vector<FloatRaster>& temperatures,
        const FloatRaster& weather_coefficient,
        Treatments<IntegerRaster, FloatRaster>& treatments,
        IntegerRaster& resistant,
        std::vector<std::tuple<int, int>>& outside_dispersers,  // out
        SpreadRate<IntegerRaster>& spread_rate,  // out
        QuarantineEscape<IntegerRaster>& quarantine,  // out
        const IntegerRaster& quarantine_areas,
        const std::vector<std::vector<int>> movements,
        const Network<RasterIndex>& network,
        std::vector<std::vector<int>>& suitable_cells,
        IntegerRaster& established_dispersers)
    {

        // removal of dispersers due to lethal temperatures
        if (config_.use_lethal_temperature && config_.lethal_schedule()[step]) {
            int lethal_step =
                simulation_step_to_action_step(config_.lethal_schedule(), step);
            simulation_.remove(
                infected,
                susceptible,
                temperatures[lethal_step],
                config_.lethal_temperature,
                suitable_cells);
        }
        // actual spread
        if (config_.spread_schedule()[step]) {
            simulation_.generate(
                dispersers,
                infected,
                config_.weather,
                weather_coefficient,
                config_.reproductive_rate,
                suitable_cells,
                established_dispersers);

            auto dispersal_kernel = kernel_factory_(config_, dispersers, network);
            auto overpopulation_kernel =
                create_overpopulation_movement_kernel(dispersers, network);

            simulation_.disperse_and_infect(
                step,
                dispersers,
                susceptible,
                exposed,
                infected,
                mortality_tracker.back(),
                total_populations,
                total_exposed,
                outside_dispersers,
                config_.weather,
                weather_coefficient,
                dispersal_kernel,
                suitable_cells,
                established_dispersers,
                config_.establishment_probability);
            if (config_.use_overpopulation_movements) {
                simulation_.move_overpopulated_pests(
                    susceptible,
                    infected,
                    total_hosts,
                    outside_dispersers,
                    overpopulation_kernel,
                    suitable_cells,
                    config_.overpopulation_percentage,
                    config_.leaving_percentage);
            }
            if (config_.use_movements) {
                // to do fix movements to use entire mortality tracker
                last_index = simulation_.movement(
                    infected,
                    susceptible,
                    mortality_tracker,
                    exposed,
                    resistant,
                    total_hosts,
                    total_exposed,
                    step,
                    last_index,
                    movements,
                    config_.movement_schedule,
                    suitable_cells);
            }
        }
        // treatments
        if (config_.use_treatments) {
            bool managed = treatments.manage(
                step,
                infected,
                exposed,
                susceptible,
                resistant,
                total_hosts,
                suitable_cells);
            if (managed && config_.use_mortality) {
                // treatments apply to all mortality tracker cohorts
                for (auto& raster : mortality_tracker) {
                    treatments.manage_mortality(step, raster, suitable_cells);
                }
            }
        }
        if (config_.use_mortality && config_.mortality_schedule()[step]) {
            // expectation is that mortality tracker is of length (1/mortality_rate
            // + mortality_time_lag).
            // TODO: died.zero(); should be done by the caller if needed, document!
            simulation_.mortality(
                infected,
                total_hosts,
                config_.mortality_rate,
                config_.mortality_time_lag,
                died,
                mortality_tracker,
                suitable_cells);
        }
        // compute spread rate
        if (config_.use_spreadrates && config_.spread_rate_schedule()[step]) {
            unsigned rates_step =
                simulation_step_to_action_step(config_.spread_rate_schedule(), step);
            spread_rate.compute_step_spread_rate(infected, rates_step, suitable_cells);
        }
        // compute quarantine escape
        if (config_.use_quarantine && config_.quarantine_schedule()[step]) {
            unsigned action_step =
                simulation_step_to_action_step(config_.quarantine_schedule(), step);
            quarantine.infection_escape_quarantine(
                infected, quarantine_areas, action_step, suitable_cells);
        }
    }
};

}  // namespace pops

#endif  // POPS_MODEL_HPP
