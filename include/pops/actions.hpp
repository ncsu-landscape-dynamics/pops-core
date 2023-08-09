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

template<
    typename Hosts,
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename Generator>
class SpreadAction
{
public:
    template<typename DispersalKernel>
    void action(
        const IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        DispersalKernel& dispersal_kernel,
        Hosts& host_pool,
        Generator& generator)
    {
        this->generate(dispersers, established_dispersers, host_pool, generator);
        this->disperse(
            dispersers,
            established_dispersers,
            outside_dispersers,
            dispersal_kernel,
            host_pool,
            generator);
    }

    /** Generates dispersers based on infected
     *
     * @param[out] dispersers  (existing values are ignored)
     * @param infected Currently infected hosts
     * @param weather Whether to use the weather coefficient
     * @param reproductive_rate reproductive rate (used unmodified when weather
     *        coefficient is not used)
     */
    void generate(
        IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        Hosts& host_pool,
        Generator& generator)
    {
        for (auto indices : host_pool.suitable_cells()) {
            int i = indices[0];
            int j = indices[1];
            if (host_pool.infected_at(i, j) > 0) {
                int dispersers_from_cell =
                    host_pool.dispersers_from(i, j, generator.disperser_generation());
                if (soil_pool_) {
                    // From all the generated dispersers, some go to the soil in the
                    // same cell and don't participate in the kernel-driven dispersal.
                    auto dispersers_to_soil =
                        std::round(to_soil_percentage_ * dispersers_from_cell);
                    soil_pool_->dispersers_to(dispersers_to_soil, i, j, generator);
                    dispersers_from_cell -= dispersers_to_soil;
                }
                dispersers(i, j) = dispersers_from_cell;
                established_dispersers(i, j) = dispersers_from_cell;
            }
            else {
                dispersers(i, j) = 0;
                established_dispersers(i, j) = 0;
            }
        }
    }

    /** Creates dispersal locations for the dispersing individuals
     *
     * Depending on what data is provided as the *exposed_or_infected*
     * parameter, this function can be part of the S to E step or the
     * S to I step.
     *
     * Typically, the generate() function is called beforehand to
     * create dispersers. In SEI model, the infect_exposed() function is
     * typically called afterwards.
     *
     * DispersalKernel is callable object or function with one parameter
     * which is the random number engine (generator). The return value
     * is row and column in the raster (or outside of it). The current
     * position is passed as parameters. The return value is in the
     * form of a tuple with row and column so that std::tie() is usable
     * on the result, i.e. function returning
     * `std::make_tuple(row, column)` fulfills this requirement.
     *
     * The *total_populations* can be total number of hosts in the basic case
     * or it can be the total size of population of all relevant species
     * both host and non-host if dilution effect should be applied.
     *
     * If establishment stochasticity is disabled,
     * *establishment_probability* is used to decide whether or not
     * a disperser is established in a cell. Value 1 means that all
     * dispersers will establish and value 0 means that no dispersers
     * will establish.
     *
     * @param[in] dispersers Dispersing individuals ready to be dispersed
     * @param[in,out] susceptible Susceptible hosts
     * @param[in,out] exposed_or_infected Exposed or infected hosts
     * @param[in,out] mortality_tracker Newly infected hosts (if applicable)
     * @param[in, out] total_exposed Total exposed in all exposed cohorts
     * @param[in] total_populations All host and non-host individuals in the area
     * @param[in,out] outside_dispersers Dispersers escaping the raster
     * @param weather Whether or not weather coefficients should be used
     * @param dispersal_kernel Dispersal kernel to move dispersers
     * @param establishment_probability Probability of establishment with no
     *        stochasticity
     *
     * @note If the parameters or their default values don't correspond
     * with the disperse_and_infect() function, it is a bug.
     */
    template<typename DispersalKernel>
    void disperse(
        const IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        DispersalKernel& dispersal_kernel,
        Hosts& host_pool,
        Generator& generator)
    {
        // The interaction does not happen over the member variables yet, use empty
        // variables. This requires SI/SEI to be fully resolved in host and not in
        // disperse_and_infect.
        int row;
        int col;
        for (auto indices : host_pool.suitable_cells()) {
            int i = indices[0];
            int j = indices[1];
            if (dispersers(i, j) > 0) {
                for (int k = 0; k < dispersers(i, j); k++) {
                    std::tie(row, col) = dispersal_kernel(generator, i, j);
                    // if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
                    if (host_pool.is_outside(row, col)) {
                        // export dispersers dispersed outside of modeled area
                        outside_dispersers.emplace_back(std::make_tuple(row, col));
                        established_dispersers(i, j) -= 1;
                        continue;
                    }
                    // Put a disperser to the host pool.
                    auto dispersed =
                        host_pool.disperser_to(row, col, generator.establishment());
                    if (!dispersed) {
                        established_dispersers(i, j) -= 1;
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
     * a percentage (0-1 ratio) of dispersers which will be send to the soil (and may
     * establish or not depending on the soil pool object).
     *
     * Soil pool is optional and implemented in more general (but experimental) way.
     * This function needs to be called separately some time after the object is created
     * to active the soil part of the simulation. This avoids the need for many
     * constructors or many optional parameters which need default values.
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
     * Optional soil pool
     */
    std::shared_ptr<SoilPool<IntegerRaster, FloatRaster, RasterIndex, Generator>>
        soil_pool_{nullptr};
    /**
     * Percentage (0-1 ratio) of disperers to be send to soil
     */
    double to_soil_percentage_{0};
};

template<typename Hosts, typename IntegerRaster, typename FloatRaster>
class SurvivalRateAction
{
public:
    SurvivalRateAction(const FloatRaster& survival_rate) : survival_rate_(survival_rate)
    {}
    /** Removes percentage of exposed and infected
     *
     * @param infected Currently infected hosts
     * @param susceptible Currently susceptible hosts
     * @param mortality_tracker_vector Hosts that are infected at a specific time
     * step
     * @param exposed Exposed hosts per cohort
     * @param total_exposed Total exposed in all exposed cohorts
     * @param survival_rate Raster between 0 and 1 representing pest survival rate
     */
    template<typename Generator>
    void action(Hosts& hosts, Generator& generator)
    {
        for (auto indices : hosts.suitable_cells()) {
            int i = indices[0];
            int j = indices[1];
            if (survival_rate_(i, j) < 1) {
                // remove percentage of infestation/infection in the infected class
                auto infected = hosts.infected_at(i, j);
                int removed_infected =
                    infected - std::lround(infected * survival_rate_(i, j));
                hosts.remove_infected_at(
                    i, j, removed_infected, generator.survival_rate());
                // remove the same percentage for total exposed and remove randomly from
                // each cohort
                auto exposed = hosts.exposed_at(i, j);
                int total_removed_exposed =
                    exposed - std::lround(exposed * survival_rate_(i, j));
                hosts.remove_exposed_at(
                    i, j, total_removed_exposed, generator.survival_rate());
            }
        }
    }

private:
    const FloatRaster& survival_rate_;
};

template<
    typename Hosts,
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename Generator>
class RemoveByTemperature
{
public:
    RemoveByTemperature(
        const EnvironmentInterface<IntegerRaster, FloatRaster, RasterIndex, Generator>&
            environment,
        double lethal_temperature)
        : environment_(environment), lethal_temperature_(lethal_temperature)
    {}
    /** removes infected based on min or max temperature tolerance
     *
     * @param infected Currently infected hosts
     * @param susceptible Currently susceptible hosts
     * @param temperature Spatially explicit temperature
     * @param lethal_temperature temperature at which lethal conditions occur
     */
    void action(
        Hosts& hosts,

        Generator& generator)
    {
        for (auto indices : hosts.suitable_cells()) {
            int i = indices[0];
            int j = indices[1];
            if (environment_.temperature_at(i, j) < lethal_temperature_) {
                auto count = hosts.infected_at(i, j);
                // TODO: Not sure what generator this should use.
                hosts.remove_infected_at(i, j, count, generator.weather());
                // now this includes also mortality, but it does not include exposed
            }
        }
    }

private:
    const EnvironmentInterface<IntegerRaster, FloatRaster, RasterIndex, Generator>&
        environment_;
    const double lethal_temperature_;
};

template<
    typename Hosts,
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex>
class MoveOverpopulatedPests
{
public:
    MoveOverpopulatedPests(
        double overpopulation_percentage,
        double leaving_percentage,
        RasterIndex rows,
        RasterIndex cols)
        : overpopulation_percentage_(overpopulation_percentage),
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
     * @param[in,out] susceptible Susceptible hosts
     * @param[in,out] infected Infected hosts
     * @param total_hosts All host individuals in the area. Is equal to
     * infected + exposed + susceptible in the cell.
     * @param[in,out] outside_dispersers Dispersers escaping the rasters
     * @param dispersal_kernel Dispersal kernel to move dispersers (pests)
     * @param overpopulation_percentage Percentage of occupied hosts when the cell is
     *        considered to be overpopulated
     * @param leaving_percentage Percentage pests leaving an overpopulated cell
     *
     * @note Exposed hosts do not count towards total number of pest,
     *       i.e., *total_host* is assumed to be S + E in SEI model.
     * @note Mortality is not supported by this function, i.e., the mortality rasters
     *       are not modified while the infected are.
     */
    template<typename DispersalKernel, typename Generator>
    void action(
        Hosts& hosts,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        DispersalKernel& dispersal_kernel,

        Generator& generator)
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
                std::tie(row, col) = dispersal_kernel(generator, i, j);
                // for leaving_percentage == 0.5
                // 2 infected -> 1 leaving
                // 3 infected -> 1 leaving
                int leaving = original_count * leaving_percentage_;
                leaving = hosts.pest_from(i, j, leaving);
                if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
                    // Collect pests dispersed outside of modeled area.
                    outside_dispersers.reserve(outside_dispersers.size() + leaving);
                    for (int pest = 0; pest < leaving; ++pest)
                        outside_dispersers.emplace_back(row, col);
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
            hosts.pests_to(move.row, move.col, move.count);
        }
    }

private:
    double overpopulation_percentage_;
    double leaving_percentage_;
    RasterIndex rows_;
    RasterIndex cols_;
};

template<
    typename Hosts,
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex>
class HostMovement
{
public:
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

    /** Moves hosts from one location to another
     *
     * @note Note that unlike the other functions, here, *total_hosts*,
     * i.e., number of hosts is required, not number of all hosts
     * and non-host individuals.
     *
     * @param infected Currently infected hosts
     * @param susceptible Currently susceptible hosts
     * @param mortality_tracker_vector Hosts that are infected at a specific time step
     * @param total_hosts All host individuals in the area. Is equal to
     *        infected + exposed + susceptible in the cell.
     * @param total_exposed Total exposed in all exposed cohorts
     * @param exposed Exposed hosts per cohort
     * @param resistant Resistant hosts
     * @param step the current step of the simulation
     * @param last_index the last index to not be used from movements
     * @param movements a vector of ints with row_from, col_from, row_to, col_to, and
     *        num_hosts
     * @param movement_schedule a vector matching movements with the step at which the
     *        movement from movements are applied
     *
     * @note Mortality and non-host individuals are not supported in movements.
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

template<typename Hosts, typename IntegerRaster, typename FloatRaster>
class Mortality
{
public:
    Mortality(double mortality_rate, int mortality_time_lag)
        : mortality_rate_(mortality_rate), mortality_time_lag_(mortality_time_lag)
    {}
    /** kills infected hosts based on mortality rate and timing. In the last year
     * of mortality tracking the first index all remaining tracked infected hosts
     * are removed. In indexes that are in the mortality_time_lag no mortality occurs.
     * In all other indexes the number of tracked individuals is multiplied by the
     * mortality rate to calculate the number of hosts that die that time step. The
     * mortality_tracker_vector has a minimum size of mortality_time_lag + 1.
     *
     * @param infected Currently infected hosts
     * @param total_hosts All hosts
     * @param mortality_rate percent of infected hosts that die each time period
     * @param mortality_time_lag time lag prior to mortality beginning
     * @param died dead hosts during time step
     * @param mortality_tracker_vector vector of matrices for tracking infected
     * host infection over time. Expectation is that mortality tracker is of
     * length (1/mortality_rate + mortality_time_lag)
     */
    void action(Hosts& hosts)
    {
        for (auto indices : hosts.suitable_cells()) {
            hosts.apply_mortality_at(
                indices[0], indices[1], mortality_rate_, mortality_time_lag_);
        }
        hosts.step_forward_mortality();
    }

private:
    const double mortality_rate_ = 0;
    const int mortality_time_lag_ = 0;
};

}  // namespace pops

#endif  // POPS_ACTIONS_HPP
