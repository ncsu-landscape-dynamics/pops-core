/*
 * PoPS model - pest or pathogen spread simulation
 *
 * Copyright (C) 2023 by the authors.
 *
 * Authors: Vaclav Petras (wenzeslaus gmail com)
 *          Chris Jones (cjones1688 gmail com)
 *          Anna Petrasova (kratochanna gmail com)
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_ACTIONS_HPP
#define POPS_ACTIONS_HPP

#include <cmath>
#include <memory>
#include <tuple>
#include <vector>
#include <random>
#include <string>
#include <stdexcept>

#include "utils.hpp"
#include "model_type.hpp"
#include "soils.hpp"

namespace pops {

/**
 * Spread of pest or pathogens
 *
 * The spread starts with generation of dispersers in the hosts.
 * This is followed by dispersion and possible establishment.
 */
template<
    typename Hosts,
    typename Pests,
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename DispersalKernel,
    typename Generator>
class SpreadAction
{
public:
    /**
     * @brief Create the object with a given dispersal kernel.
     *
     * DispersalKernel is callable object or function with one parameter
     * which is the random number engine (generator). The return value
     * is row and column in the raster (or outside of it). The current
     * position is passed as parameters. The return value is in the
     * form of a tuple with row and column so that std::tie() is usable
     * on the result, i.e. function returning
     * `std::make_tuple(row, column)` fulfills this requirement.
     *
     * @param dispersal_kernel Dispersal kernel to move the dispersers
     */
    SpreadAction(DispersalKernel& dispersal_kernel)
        : dispersal_kernel_(dispersal_kernel)
    {}

    /** Perform the generation and spread */
    void action(Hosts& host_pool, Pests& pests, Generator& generator)
    {
        this->generate(host_pool, pests, generator);
        this->disperse(host_pool, pests, generator);
    }

    /** Generates dispersers based on hosts
     *
     * @note Unlike other functions, this resets the state of dispersers in the pest
     * pool.
     */
    void generate(Hosts& host_pool, Pests& pests, Generator& generator)
    {
        for (auto indices : host_pool.suitable_cells()) {
            int i = indices[0];
            int j = indices[1];
            int dispersers_from_cell =
                host_pool.dispersers_from(i, j, generator.disperser_generation());
            if (dispersers_from_cell > 0) {
                if (soil_pool_) {
                    // From all the generated dispersers, some go to the soil in the
                    // same cell and don't participate in the kernel-driven dispersal.
                    auto dispersers_to_soil =
                        std::lround(to_soil_percentage_ * dispersers_from_cell);
                    soil_pool_->dispersers_to(dispersers_to_soil, i, j, generator);
                    dispersers_from_cell -= dispersers_to_soil;
                }
                pests.set_dispersers_at(i, j, dispersers_from_cell, 0);
            }
            else {
                pests.set_dispersers_at(i, j, 0, 0);
            }
        }
    }

    /** Moves dispersers (dispersing individuals) to dispersal locations
     *
     * The function sends dispersers previously generated by generate() to locations in
     * host pool generated by the dispersal kernel. Any SI-SEI differences are resolved
     * by the host.
     */
    void disperse(Hosts& host_pool, Pests& pests, Generator& generator)
    {
        // The interaction does not happen over the member variables yet, use empty
        // variables. This requires SI/SEI to be fully resolved in host and not in
        // disperse_and_infect.
        int row;
        int col;
        for (auto indices : host_pool.suitable_cells()) {
            int i = indices[0];
            int j = indices[1];
            if (pests.dispersers_at(i, j) > 0) {
                for (int k = 0; k < pests.dispersers_at(i, j); k++) {
                    std::tie(row, col) = dispersal_kernel_(generator, i, j);
                    if (host_pool.is_outside(row, col)) {
                        pests.add_outside_disperser_at(row, col);
                        continue;
                    }
                    // Put a disperser to the host pool.
                    auto dispersed =
                        host_pool.disperser_to(row, col, generator.establishment());
                    if (dispersed) {
                        pests.add_established_dispersers_at(i, j, 1);
                    }
                }
            }
            if (soil_pool_) {
                // Get dispersers from the soil if there are any.
                auto num_dispersers = soil_pool_->dispersers_from(i, j, generator);
                // Put each disperser to the host pool.
                for (int k = 0; k < num_dispersers; k++) {
                    host_pool.disperser_to(i, j, generator.establishment());
                }
            }
        }
    }

    /**
     * @brief Activate storage of dispersers in soil
     *
     * Calling this function activates the soils. By default, the soil pool is not used.
     * The parameters are soil pool used to store the dispersers and
     * a percentage (0-1 ratio) of dispersers which will be sent to the soil (and may
     * establish or not depending on the soil pool object).
     *
     * Soil pool is optional and nothing is done when it is not set.
     * This function needs to be called separately some time after the object is created
     * to activate the soil part of the simulation. This avoids the need for many
     * more constructors or for many optional parameters which need default values.
     *
     * @param soil_pool Soils pool object to use for storage
     * @param dispersers_percentage Percentage of dispersers moving to the soil
     */
    void activate_soils(
        std::shared_ptr<SoilPool<IntegerRaster, FloatRaster, RasterIndex, Generator>>
            soil_pool,
        double dispersers_percentage)
    {
        this->soil_pool_ = soil_pool;
        this->to_soil_percentage_ = dispersers_percentage;
    }

private:
    /**
     * Dispersal kernel
     */
    DispersalKernel& dispersal_kernel_;
    /**
     * Optional soil pool
     */
    std::shared_ptr<SoilPool<IntegerRaster, FloatRaster, RasterIndex, Generator>>
        soil_pool_{nullptr};
    /**
     * Percentage (0-1 ratio) of disperers to be send to soil
     */
    double to_soil_percentage_{0};
};

/**
 * Pest or pathogen individuals survive only with a given survival rate
 *
 * This removes a given percentage of exposed and infected in host pool.
 */
template<typename Hosts, typename IntegerRaster, typename FloatRaster>
class SurvivalRateAction
{
public:
    /**
     * @brief Create object with a given survival rate.
     *
     * Survival rate is a spatially variable float raster with values between 0 and 1.
     * From pest perspective this is the ratio of surviving individuals on hosts,
     * while from the host perspective this is the ratio of hosts to keep infected or
     * exposed.
     *
     * @param survival_rate Survival rate
     */
    SurvivalRateAction(const FloatRaster& survival_rate) : survival_rate_(survival_rate)
    {}
    /**
     * Reduce the infection based on the pest survival rate.
     *
     * Infected and total number of exposed hosts are removed directly, while mortality
     * cohorts and exposed cohorts are left to be managed by the host pool. In other
     * words, the details of how total infected and total exposed are distributed among
     * the cohorts is managed by the host pools.
     */
    template<typename Generator>
    void action(Hosts& hosts, Generator& generator)
    {
        for (auto indices : hosts.suitable_cells()) {
            int i = indices[0];
            int j = indices[1];
            if (survival_rate_(i, j) < 1) {
                hosts.remove_infection_by_ratio_at(
                    i, j, survival_rate_(i, j), generator.survival_rate());
            }
        }
    }

private:
    const FloatRaster& survival_rate_;
};

/**
 * Remove infection or infestation by temperature threshold
 */
template<
    typename Hosts,
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename Generator>
class RemoveByTemperature
{
public:
    /**
     * @brief Create an object with reference to the environment and the temperature
     * threshold
     *
     * All infection in a cell with temperature under the *lethal_temperature* value is
     * removed from hosts.
     *
     * @param environment Environment which supplies the temperature per cell
     * @param lethal_temperature Temperature at which lethal conditions occur
     */
    RemoveByTemperature(
        const EnvironmentInterface<IntegerRaster, FloatRaster, RasterIndex, Generator>&
            environment,
        double lethal_temperature)
        : environment_(environment), lethal_temperature_(lethal_temperature)
    {}
    /** Perform the removal of infection */
    void action(Hosts& hosts, Generator& generator)
    {
        for (auto indices : hosts.suitable_cells()) {
            int i = indices[0];
            int j = indices[1];
            if (environment_.temperature_at(i, j) < lethal_temperature_) {
                // now this includes also mortality, but it does not include exposed
                hosts.remove_all_infected_at(i, j, generator.lethal_temperature());
            }
        }
    }

private:
    const EnvironmentInterface<IntegerRaster, FloatRaster, RasterIndex, Generator>&
        environment_;
    const double lethal_temperature_;
};

/**
 * Move overpopulated pests in one cell to a different cell.
 *
 * Overpopulation is measured by the ratio of infected hosts to all hosts.
 */
template<
    typename Hosts,
    typename Pests,
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename DispersalKernel>
class MoveOverpopulatedPests
{
public:
    /**
     * @brief Create object with given dispersal kernel and overpopulation parameters
     *
     * @param dispersal_kernel Dispersal kernel used for pest individuals leaving a cell
     * @param overpopulation_percentage Ratio of infected to all hosts considered as
     *        overpopulation
     * @param leaving_percentage Ratio of pest individuals leaving the cell
     * @param rows Number of rows in the whole area (for tracking outside dispersers)
     * @param cols Number of columns in the whole area (for tracking outside dispersers)
     */
    MoveOverpopulatedPests(
        DispersalKernel& dispersal_kernel,
        double overpopulation_percentage,
        double leaving_percentage,
        RasterIndex rows,
        RasterIndex cols)
        : dispersal_kernel_(dispersal_kernel),
          overpopulation_percentage_(overpopulation_percentage),
          leaving_percentage_(leaving_percentage),
          rows_(rows),
          cols_(cols)
    {}
    /** Move overflowing pest population to other hosts.
     *
     * When the number of pests (pest population) is too high, part of them moves
     * to a different location. Number of infected/infested hosts is considered to be
     * the number of pests (groups of pest) in a raster cell.
     *
     * The movement happens in two stages. First, all the leaving pests are identified
     * and removed from the source cells. Second, the move to the target cells is
     * performed. This means that even if the resulting number of pests in the target
     * cell is considered too high, it is left as is and the move is performed the next
     * time this function is called.
     *
     * If the pests (pest population) cannot be accommodated in the target cell due to
     * the insufficient number of susceptible hosts, the excessive pests die.
     *
     * @note Exposed hosts do not count towards total number of pest,
     *       i.e., *total_host* is assumed to be S + E in SEI model.
     * @note Mortality and exposure are not supported by this function, i.e., the
     * mortality rasters are not modified while the infected are. This is left to be
     * managed by the host pool, but the current host pool does not support that.
     */
    template<typename Generator>
    void action(Hosts& hosts, Pests& pests, Generator& generator)
    {
        struct Move
        {
            RasterIndex row;
            RasterIndex col;
            int count;
        };
        std::vector<Move> moves;

        // Identify the moves. Remove from source cells.
        for (auto indices : hosts.suitable_cells()) {
            int i = indices[0];
            int j = indices[1];
            int original_count = hosts.infected_at(i, j);
            // No move with only one infected host (one unit).
            if (original_count <= 1)
                continue;
            // r = I / (I + S)
            // r = I / (I + S + E_1 + E_2 + ...)
            double ratio = original_count / double(hosts.total_hosts_at(i, j));
            if (ratio >= overpopulation_percentage_) {
                int row;
                int col;
                std::tie(row, col) =
                    dispersal_kernel_(generator.overpopulation(), i, j);
                // for leaving_percentage == 0.5
                // 2 infected -> 1 leaving
                // 3 infected -> 2 leaving (assuming always rounding up for .5)
                int leaving = std::lround(original_count * leaving_percentage_);
                leaving = hosts.pests_from(i, j, leaving, generator.overpopulation());
                if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
                    pests.add_outside_dispersers_at(row, col, leaving);
                    continue;
                }
                // Doing the move here would create inconsistent results as some
                // target cells would be evaluated after the moved pest arrived,
                // possibly triggering another move.
                // So, instead, we just collect them and apply later. (The pest is in
                // the air when in the moves vector.)
                moves.push_back(Move({row, col, leaving}));
            }
        }
        // Perform the moves to target cells.
        for (const auto& move : moves) {
            // Pests which won't fit are ignored (disappear). This can happen if there
            // is simply not enough S hosts to accommodate all the pests moving from the
            // source or if multiple sources end up in the same target cell and there is
            // not enough S hosts to accommodate all of them. The decision is made in
            // the host pool. Here, we ignore the return value specifying the number of
            // accepted pests.
            hosts.pests_to(move.row, move.col, move.count, generator.overpopulation());
        }
    }

private:
    DispersalKernel& dispersal_kernel_;
    double overpopulation_percentage_;
    double leaving_percentage_;
    RasterIndex rows_;
    RasterIndex cols_;
};

/**
 * Moves hosts from one location to another
 *
 * @note Mortality and non-host individuals are not supported in movements.
 */
template<
    typename Hosts,
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex>
class HostMovement
{
public:
    /**
     * @brief Creates an object with state pointing to the movements to be used
     *
     * @param step the current step of the simulation
     * @param last_index the last index to not be used from movements
     * @param movements a vector of ints with row_from, col_from, row_to, col_to, and
     *        num_hosts
     * @param movement_schedule a vector matching movements with the step at which the
     *        movement from movements are applied
     *
     * @note There is a mix of signed and unsigned ints for step and movement schedule.
     */
    HostMovement(
        unsigned step,
        unsigned last_index,
        const std::vector<std::vector<int>>& movements,
        const std::vector<unsigned>& movement_schedule)
        : step_(step),
          last_index_(last_index),
          movements_(movements),
          movement_schedule_(movement_schedule)
    {}

    /** Perform the movement
     *
     * @note This one returns a value, so to create a consistent interface for all
     * actions, the tracking would have to happen in the action which would be better
     * design anyway.
     */
    template<typename Generator>
    unsigned action(Hosts& hosts, Generator& generator)
    {
        for (unsigned i = last_index_; i < movements_.size(); i++) {
            const auto& moved = movements_[i];
            unsigned move_schedule = movement_schedule_[i];
            if (move_schedule != step_) {
                return i;
            }
            int row_from = moved[0];
            int col_from = moved[1];
            int row_to = moved[2];
            int col_to = moved[3];
            hosts.move_hosts_from_to(
                row_from, col_from, row_to, col_to, moved[4], generator.movement());
        }
        return movements_.size();
    }

private:
    const unsigned step_;
    const unsigned last_index_;
    const std::vector<std::vector<int>>& movements_;
    const std::vector<unsigned>& movement_schedule_;
};

/**
 * Mortality of the hosts
 *
 * Kills infected hosts based on mortality rate and timing. The host pool implementation
 * of mortality is used to perform the actual process while this class provides
 * parameters and takes care of the two steps in the mortality process (dying and
 * aging).
 *
 * The *mortality_tracker_vector* used by hosts must fit with the *mortality_time_lag*
 * here. It needs to have a minimum size of mortality_time_lag + 1. See
 * HostPool::apply_mortality_at() and HostPool::step_forward_mortality().
 */
template<typename Hosts, typename IntegerRaster, typename FloatRaster>
class Mortality
{
public:
    /**
     * @brief Create object which will let host decide its mortality parameters.
     */
    Mortality() : action_mortality_(false) {}
    /**
     * @brief Create object with fixed mortality rate and time lag.
     *
     * @param mortality_rate Percent of infected hosts that die each time period
     * @param mortality_time_lag Time lag prior to mortality beginning
     */
    Mortality(double mortality_rate, int mortality_time_lag)
        : mortality_rate_(mortality_rate),
          mortality_time_lag_(mortality_time_lag),
          action_mortality_(true)
    {}
    /**
     * Perform the action by applying mortality and moving the mortality tracker
     * forward.
     */
    void action(Hosts& hosts)
    {
        for (auto indices : hosts.suitable_cells()) {
            if (static_cast<bool>(action_mortality_)) {
                hosts.apply_mortality_at(
                    indices[0], indices[1], mortality_rate_, mortality_time_lag_);
            }
            else {
                hosts.apply_mortality_at(indices[0], indices[1]);
            }
        }
        hosts.step_forward_mortality();
    }

private:
    const double mortality_rate_ = 0;
    const int mortality_time_lag_ = 0;
    const double action_mortality_ = false;
};

}  // namespace pops

#endif  // POPS_ACTIONS_HPP
