/*
 * PoPS model - pest or pathogen spread simulation
 *
 * Copyright (C) 2015-2022 by the authors.
 *
 * Authors: Zexi Chen (zchen22 ncsu edu)
 *          Vaclav Petras (wenzeslaus gmail com)
 *          Anna Petrasova (kratochanna gmail com)
 *          Chris Jones (cjones1688 gmail com)
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_SIMULATION_HPP
#define POPS_SIMULATION_HPP

#include <cmath>
#include <memory>
#include <tuple>
#include <vector>
#include <random>
#include <string>
#include <stdexcept>

#include "utils.hpp"
#include "soils.hpp"
#include "model_type.hpp"
#include "host_pool.hpp"

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
        const std::vector<std::vector<int>>& suitable_cells,
        Generator& generator)
    {
        this->generate(
            dispersers, established_dispersers, host_pool, suitable_cells, generator);
        this->disperse(
            dispersers,
            established_dispersers,
            outside_dispersers,
            dispersal_kernel,
            host_pool,
            suitable_cells,
            generator);
    }

    /** Generates dispersers based on infected
     *
     * @param[out] dispersers  (existing values are ignored)
     * @param infected Currently infected hosts
     * @param weather Whether to use the weather coefficient
     * @param reproductive_rate reproductive rate (used unmodified when weather
     *        coefficient is not used)
     * @param[in] suitable_cells List of indices of cells with hosts
     */
    void generate(
        IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        Hosts& host_pool,
        const std::vector<std::vector<int>>& suitable_cells,
        Generator& generator)
    {
        for (auto indices : suitable_cells) {
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
     * @param[in] suitable_cells List of indices of cells with hosts
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
        const std::vector<std::vector<int>>& suitable_cells,
        Generator& generator)
    {
        // The interaction does not happen over the member variables yet, use empty
        // variables. This requires SI/SEI to be fully resolved in host and not in
        // disperse_and_infect.
        int row;
        int col;
        for (auto indices : suitable_cells) {
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
     * @param suitable_cells used to run model only where host are known to occur
     */
    template<typename Generator>
    void action(
        Hosts& hosts,
        const std::vector<std::vector<int>>& suitable_cells,
        Generator& generator)
    {
        for (auto indices : suitable_cells) {
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

template<typename Hosts, typename IntegerRaster, typename FloatRaster>
class RemoveByTemperature
{
public:
    /** removes infected based on min or max temperature tolerance
     *
     * @param infected Currently infected hosts
     * @param susceptible Currently susceptible hosts
     * @param temperature Spatially explicit temperature
     * @param lethal_temperature temperature at which lethal conditions occur
     * @param suitable_cells used to run model only where host are known to occur
     */
    template<typename Generator>
    void action(
        Hosts& hosts,
        const FloatRaster& temperature,
        double lethal_temperature,
        const std::vector<std::vector<int>>& suitable_cells,
        Generator& generator)
    {
        for (auto indices : suitable_cells) {
            int i = indices[0];
            int j = indices[1];
            if (temperature(i, j) < lethal_temperature) {
                auto count = hosts.infected_at(i, j);
                // TODO: Not sure what generator this should use.
                hosts.remove_infected_at(i, j, count, generator.weather());
                // now this includes also mortality, but it does not include exposed
            }
        }
    }
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
     * @param[in] suitable_cells List of indices of cells with hosts
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
        const std::vector<std::vector<int>>& suitable_cells,
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
        for (auto indices : suitable_cells) {
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
     * @param suitable_cells List of indices of cells with hosts
     *
     * @note Mortality and non-host individuals are not supported in movements.
     */
    template<typename Generator>
    unsigned movement(
        Hosts& hosts,
        unsigned step,
        unsigned last_index,
        const std::vector<std::vector<int>>& movements,
        std::vector<unsigned> movement_schedule,
        Generator& generator)
    {
        for (unsigned i = last_index; i < movements.size(); i++) {
            auto moved = movements[i];
            unsigned move_schedule = movement_schedule[i];
            if (move_schedule != step) {
                return i;
            }
            int row_from = moved[0];
            int col_from = moved[1];
            int row_to = moved[2];
            int col_to = moved[3];
            hosts.move_hosts_from_to(
                row_from, col_from, row_to, col_to, moved[4], generator.movement());
        }
        return movements.size();
    }
};

template<typename Hosts, typename IntegerRaster, typename FloatRaster>
class Mortality
{
public:
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
     * @param suitable_cells used to run model only where host are known to occur
     */
    void action(
        Hosts& hosts,
        double mortality_rate,
        int mortality_time_lag,
        const std::vector<std::vector<int>>& suitable_cells)
    {
        for (auto indices : suitable_cells) {
            hosts.apply_mortality_at(
                indices[0], indices[1], mortality_rate, mortality_time_lag);
        }
        hosts.step_forward_mortality();
    }
};

/*! The main class to control the spread simulation.
 *
 * The Simulation class handles the mechanics of the model, but the
 * timing of the events or steps should be handled outside of this
 * class unless noted otherwise. The notable exceptions are exposed
 * hosts in the SEI model type and mortality.
 *
 * The template parameters IntegerRaster and FloatRaster are raster
 * image or matrix types. Any 2D numerical array should work as long as
 * it uses function call operator to access the values, i.e. it provides
 * indexing for reading and writing values using `()`. In other words,
 * the operations such as the two following ones should be possible:
 *
 * ```
 * a(i, j) = 1;
 * a(i, j) == 1;
 * ```
 *
 * The PoPS library offers a Raster template class to fill this role,
 * but other classes can be used as well.
 *
 * Template parameter RasterIndex is type used for maximum indices of
 * the used rasters and should be the same as what the actual raster
 * types are using. However, at the same time, comparison with signed
 * type are performed and a signed type might be required in the future.
 * A default is provided, but it can be changed in the future.
 */
template<
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex = int,
    typename Generator = DefaultSingleGeneratorProvider>
class Simulation
{
private:
    RasterIndex rows_;
    RasterIndex cols_;
    bool dispersers_stochasticity_;
    bool establishment_stochasticity_;
    bool movement_stochasticity_;
    ModelType model_type_;
    unsigned latency_period_;
    /// Non-owning pointer to environment for weather
    // can be const in simulation, except for need to add host now
    Environment<IntegerRaster, FloatRaster, RasterIndex, Generator>* environment_{
        nullptr};
    /**
     * Optional soil pool
     */
    std::shared_ptr<SoilPool<IntegerRaster, FloatRaster, RasterIndex, Generator>>
        soil_pool_{nullptr};
    /**
     * Percentage (0-1 ratio) of disperers to be send to soil
     */
    double to_soil_percentage_{0};

public:
    // Host pool has the provider from model, but in test, it gets plain engine.
    using StandardHostPool =
        HostPool<IntegerRaster, FloatRaster, RasterIndex, Generator>;

    /** Creates simulation object and seeds the internal random number generator.
     *
     * The same random number generator is used throughout the simulation
     * and is seeded once at the beginning.
     *
     * The number or rows and columns needs to be the same as the size
     * of rasters used with the Simulation object
     * (potentially, it can be also smaller).
     *
     * @param model_type Type of the model (SI or SEI)
     * @param latency_period Length of the latency period in steps
     * @param random_seed Number to seed the random number generator
     * @param rows Number of rows
     * @param cols Number of columns
     * @param dispersers_stochasticity Enable stochasticity in generating of dispersers
     * @param establishment_stochasticity Enable stochasticity in establishment step
     * @param movement_stochasticity Enable stochasticity in movement of hosts
     */
    Simulation(
        RasterIndex rows,
        RasterIndex cols,
        ModelType model_type = ModelType::SusceptibleInfected,
        unsigned latency_period = 0,
        bool dispersers_stochasticity = true,
        bool establishment_stochasticity = true,
        bool movement_stochasticity = true)
        : rows_(rows),
          cols_(cols),
          dispersers_stochasticity_(dispersers_stochasticity),
          establishment_stochasticity_(establishment_stochasticity),
          movement_stochasticity_(movement_stochasticity),
          model_type_(model_type),
          latency_period_(latency_period)
    {}

    Simulation() = delete;

    /**
     * @brief Set environment used for weather to provided environment
     * @param environment Pointer to an existing environment
     *
     * The simulation object does not take ownership of the environment.
     */
    // parameter and attribute should be const
    void set_environment(
        Environment<IntegerRaster, FloatRaster, RasterIndex, Generator>* environment)
    {
        this->environment_ = environment;
    }

    /**
     * @brief Get environment used in the simulation
     *
     * @param allow_empty if true, empty (non-functional) environment is returned
     * @return Const pointer to the environment
     * @throw std::logic_error when environment is not set
     */
    const Environment<IntegerRaster, FloatRaster, RasterIndex, Generator>*
    environment(bool allow_empty = false)
    {
        static Environment<IntegerRaster, FloatRaster, RasterIndex, Generator> empty;
        if (!this->environment_) {
            if (allow_empty)
                return &empty;
            throw std::logic_error("Environment used in Simulation, but not provided");
        }
        return this->environment_;
    }

    /** removes infected based on min or max temperature tolerance
     *
     * @param infected Currently infected hosts
     * @param susceptible Currently susceptible hosts
     * @param temperature Spatially explicit temperature
     * @param lethal_temperature temperature at which lethal conditions occur
     * @param suitable_cells used to run model only where host are known to occur
     */
    template<typename GeneratorProvider>
    void remove(
        IntegerRaster& infected,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& total_exposed,
        std::vector<IntegerRaster>& mortality_tracker_vector,
        const FloatRaster& temperature,
        double lethal_temperature,
        std::vector<std::vector<int>>& suitable_cells,
        GeneratorProvider& generator)
    {
        IntegerRaster empty;
        StandardHostPool hosts(
            model_type_,
            susceptible,
            exposed,
            0,
            infected,
            total_exposed,
            empty,
            mortality_tracker_vector,
            empty,
            empty,
            *environment(true),
            false,
            0,
            false,
            0,
            0,
            0,
            suitable_cells);
        RemoveByTemperature<StandardHostPool, IntegerRaster, FloatRaster> remove;
        remove.action(
            hosts, temperature, lethal_temperature, suitable_cells, generator);
    }

    /** Removes percentage of exposed and infected
     *
     * @param infected Currently infected hosts
     * @param susceptible Currently susceptible hosts
     * @param mortality_tracker_vector Hosts that are infected at a specific time step
     * @param exposed Exposed hosts per cohort
     * @param total_exposed Total exposed in all exposed cohorts
     * @param survival_rate Raster between 0 and 1 representing pest survival rate
     * @param suitable_cells used to run model only where host are known to occur
     */
    template<typename GeneratorProvider>
    void remove_percentage(
        IntegerRaster& infected,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& mortality_tracker_vector,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& total_exposed,
        const FloatRaster& survival_rate,
        std::vector<std::vector<int>>& suitable_cells,
        GeneratorProvider& generator)
    {
        IntegerRaster empty;
        StandardHostPool hosts(
            model_type_,
            susceptible,
            exposed,
            0,
            infected,
            total_exposed,
            empty,
            mortality_tracker_vector,
            empty,
            empty,
            *environment(true),
            false,
            0,
            false,
            0,
            0,
            0,
            suitable_cells);
        SurvivalRateAction<StandardHostPool, IntegerRaster, FloatRaster> survival(
            survival_rate);
        survival.action(hosts, suitable_cells, generator);
    }

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
     * @param suitable_cells used to run model only where host are known to occur
     */
    void mortality(
        IntegerRaster& infected,
        IntegerRaster& total_hosts,
        double mortality_rate,
        int mortality_time_lag,
        IntegerRaster& died,
        std::vector<IntegerRaster>& mortality_tracker_vector,
        std::vector<std::vector<int>>& suitable_cells)
    {
        IntegerRaster empty;
        std::vector<IntegerRaster> empty_vector;
        StandardHostPool hosts{
            model_type_,
            empty,
            empty_vector,
            0,
            infected,
            empty,
            empty,
            mortality_tracker_vector,
            died,
            total_hosts,
            *environment(true),
            false,
            0,
            false,
            0,
            0,
            0,
            suitable_cells};
        Mortality<StandardHostPool, IntegerRaster, FloatRaster> mortality;
        mortality.action(hosts, mortality_rate, mortality_time_lag, suitable_cells);
    }

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
     * @param suitable_cells List of indices of cells with hosts
     *
     * @note Mortality and non-host individuals are not supported in movements.
     */
    unsigned movement(
        IntegerRaster& infected,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& mortality_tracker_vector,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& resistant,
        IntegerRaster& total_hosts,
        IntegerRaster& total_exposed,
        unsigned step,
        unsigned last_index,
        const std::vector<std::vector<int>>& movements,
        std::vector<unsigned> movement_schedule,
        std::vector<std::vector<int>>& suitable_cells,
        Generator& generator)
    {
        HostMovement<StandardHostPool, IntegerRaster, FloatRaster, RasterIndex>
            host_movement{};
        IntegerRaster empty;
        StandardHostPool hosts{
            model_type_,
            susceptible,
            exposed,
            0,
            infected,
            total_exposed,
            resistant,
            mortality_tracker_vector,
            empty,
            total_hosts,
            *environment(true),
            false,
            0,
            false,
            0,
            0,
            0,
            suitable_cells};
        return host_movement.movement(
            hosts, step, last_index, movements, movement_schedule, generator);
    }

    /** Generates dispersers based on infected
     *
     * @param[out] dispersers  (existing values are ignored)
     * @param infected Currently infected hosts
     * @param weather Whether to use the weather coefficient
     * @param reproductive_rate reproductive rate (used unmodified when weather
     *        coefficient is not used)
     * @param[in] suitable_cells List of indices of cells with hosts
     */
    void generate(
        IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        const IntegerRaster& infected,
        bool weather,
        double reproductive_rate,
        const std::vector<std::vector<int>>& suitable_cells,
        Generator& generator)
    {
        IntegerRaster empty;
        std::vector<IntegerRaster> empty_vector;
        StandardHostPool host_pool{
            model_type_,
            empty,
            empty_vector,
            0,
            const_cast<IntegerRaster&>(infected),
            empty,
            empty,
            empty_vector,
            empty,
            empty,
            *environment(!weather),
            dispersers_stochasticity_,
            reproductive_rate,
            false,
            0,
            0,
            0,
            const_cast<std::vector<std::vector<int>>&>(suitable_cells)};
        SpreadAction<
            StandardHostPool,
            IntegerRaster,
            FloatRaster,
            RasterIndex,
            Generator>
            spread_action;
        spread_action.activate_soils(soil_pool_, to_soil_percentage_);
        spread_action.generate(
            dispersers, established_dispersers, host_pool, suitable_cells, generator);
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
     * @param[in] suitable_cells List of indices of cells with hosts
     *
     * @note If the parameters or their default values don't correspond
     * with the disperse_and_infect() function, it is a bug.
     */
    template<typename DispersalKernel>
    void disperse(
        const IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& infected,
        std::vector<IntegerRaster>& mortality_tracker,
        const IntegerRaster& total_populations,
        IntegerRaster& total_exposed,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        bool weather,
        DispersalKernel& dispersal_kernel,
        std::vector<std::vector<int>>& suitable_cells,
        double establishment_probability,
        Generator& generator)
    {
        // The interaction does not happen over the member variables yet, use empty
        // variables. This requires SI/SEI to be fully resolved in host and not in
        // disperse_and_infect.
        IntegerRaster empty;
        StandardHostPool host_pool{
            model_type_,
            susceptible,
            exposed,
            0,
            infected,
            total_exposed,
            empty,
            mortality_tracker,
            empty,
            empty,
            *environment(!weather),
            false,
            0,
            establishment_stochasticity_,
            establishment_probability,
            rows_,
            cols_,
            suitable_cells};
        // This would be part of the main initialization process.
        if (environment_) {
            environment_->set_total_population(&total_populations);
        }

        SpreadAction<
            StandardHostPool,
            IntegerRaster,
            FloatRaster,
            RasterIndex,
            Generator>
            spread_action;
        spread_action.activate_soils(soil_pool_, to_soil_percentage_);
        spread_action.disperse(
            dispersers,
            established_dispersers,
            outside_dispersers,
            dispersal_kernel,
            host_pool,
            suitable_cells,
            generator);
    }

    // For backwards compatibility for tests (without exposed and mortality)
    template<typename DispersalKernel>
    void disperse(
        const IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        IntegerRaster& susceptible,
        IntegerRaster& infected,
        IntegerRaster& mortality_tracker,
        const IntegerRaster& total_populations,
        IntegerRaster& total_exposed,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        bool weather,
        DispersalKernel& dispersal_kernel,
        std::vector<std::vector<int>>& suitable_cells,
        double establishment_probability = 0.5)
    {
        std::vector<IntegerRaster> tmp;
        tmp.push_back(mortality_tracker);
        std::vector<IntegerRaster> empty_vector;
        disperse(
            dispersers,
            established_dispersers,
            susceptible,
            empty_vector,
            infected,
            tmp,  // mortality
            total_populations,
            total_exposed,
            outside_dispersers,
            weather,
            dispersal_kernel,
            suitable_cells,
            establishment_probability);
        mortality_tracker = tmp.back();
    }

    // For backwards compatibility for tests (without exposed and mortality)
    template<typename DispersalKernel>
    void disperse(
        const IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        IntegerRaster& susceptible,
        IntegerRaster& infected,
        IntegerRaster& mortality_tracker,
        const IntegerRaster& total_populations,
        IntegerRaster& total_exposed,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        bool weather,
        DispersalKernel& dispersal_kernel,
        std::vector<std::vector<int>>& suitable_cells,
        double establishment_probability,
        Generator& generator)
    {
        std::vector<IntegerRaster> tmp;
        tmp.push_back(mortality_tracker);
        std::vector<IntegerRaster> empty_vector;
        disperse(
            dispersers,
            established_dispersers,
            susceptible,
            empty_vector,
            infected,
            tmp,  // mortality
            total_populations,
            total_exposed,
            outside_dispersers,
            weather,
            dispersal_kernel,
            suitable_cells,
            establishment_probability,
            generator);
        mortality_tracker = tmp.back();
    }

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
     * @param[in] suitable_cells List of indices of cells with hosts
     * @param overpopulation_percentage Percentage of occupied hosts when the cell is
     *        considered to be overpopulated
     * @param leaving_percentage Percentage pests leaving an overpopulated cell
     *
     * @note Exposed hosts do not count towards total number of pest,
     *       i.e., *total_host* is assumed to be S + E in SEI model.
     * @note Mortality is not supported by this function, i.e., the mortality rasters
     *       are not modified while the infected are.
     */
    template<typename DispersalKernel>
    void move_overpopulated_pests(
        IntegerRaster& susceptible,
        IntegerRaster& infected,
        const IntegerRaster& total_hosts,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        DispersalKernel& dispersal_kernel,
        std::vector<std::vector<int>>& suitable_cells,
        double overpopulation_percentage,
        double leaving_percentage,
        Generator& generator)
    {
        UNUSED(total_hosts);  // Total hosts is computed now.
        IntegerRaster empty;
        std::vector<IntegerRaster> empty_vector;
        StandardHostPool hosts{
            model_type_,
            susceptible,
            empty_vector,
            0,
            infected,
            empty,
            empty,
            empty_vector,
            empty,
            empty,
            *environment(true),
            false,
            0,
            false,
            0,
            0,
            0,
            suitable_cells};
        MoveOverpopulatedPests<
            StandardHostPool,
            IntegerRaster,
            FloatRaster,
            RasterIndex>
            move_pest{overpopulation_percentage, leaving_percentage, rows_, cols_};
        move_pest.action(
            hosts, outside_dispersers, dispersal_kernel, suitable_cells, generator);
    }

    /** Disperse, expose, and infect based on dispersers
     *
     * This function wraps disperse() and infect_exposed() for use in SI
     * and SEI models.
     *
     * In case of SEI model, before calling this function, last item in
     * the exposed vector needs to be ready to be used for exposure,
     * i.e., typically, it should be empty in the sense that there are
     * no hosts in the raster. This is normally taken care of by a
     * previous call to this function. The initial state of the exposed
     * vector should be such that size is latency period in steps plus 1
     * and each raster is empty, i.e., does not contain any hosts
     * (all values set to zero).
     *
     * See the infect_exposed() function for the details about exposed
     * vector, its size, and its items.
     *
     * See disperse() and infect_exposed() for a detailed list of
     * parameters and behavior. The disperse() parameter documentation
     * can be applied as is except that disperse() function's parameter
     * *exposed_or_infested* is expected to change based on the context
     * while this function's parameter *infected* is always the infected
     * individuals. Besides parameters from disperse(), this function
     * has parameter *exposed* which is the same as the one from the
     * infect_exposed() function.
     */
    template<typename DispersalKernel>
    void disperse_and_infect(
        unsigned step,
        const IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& infected,
        std::vector<IntegerRaster>& mortality_tracker,
        const IntegerRaster& total_populations,
        IntegerRaster& total_exposed,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        bool weather,
        DispersalKernel& dispersal_kernel,
        std::vector<std::vector<int>>& suitable_cells,
        double establishment_probability,
        Generator& generator)
    {
        this->disperse(
            dispersers,
            established_dispersers,
            susceptible,
            exposed,
            infected,
            mortality_tracker,
            total_populations,
            total_exposed,
            outside_dispersers,
            weather,
            dispersal_kernel,
            suitable_cells,
            establishment_probability,
            generator);
        IntegerRaster empty;
        StandardHostPool host_pool{
            model_type_,
            susceptible,
            exposed,
            latency_period_,
            infected,
            total_exposed,
            empty,
            mortality_tracker,
            empty,
            empty,
            *environment(!weather),
            false,
            0,
            establishment_stochasticity_,
            establishment_probability,
            rows_,
            cols_,
            suitable_cells};
        host_pool.step_forward(step);
    }

    template<typename DispersalKernel>
    void disperse_and_infect(
        unsigned step,
        const IntegerRaster& dispersers,
        IntegerRaster& established_dispersers,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& infected,
        IntegerRaster& mortality_tracker,
        const IntegerRaster& total_populations,
        IntegerRaster& total_exposed,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        bool weather,
        DispersalKernel& dispersal_kernel,
        std::vector<std::vector<int>>& suitable_cells,
        double establishment_probability,
        Generator& generator)
    {
        std::vector<IntegerRaster> tmp;
        tmp.push_back(mortality_tracker);
        disperse_and_infect(
            step,
            dispersers,
            established_dispersers,
            susceptible,
            exposed,
            infected,
            tmp,
            total_populations,
            total_exposed,
            outside_dispersers,
            weather,
            dispersal_kernel,
            suitable_cells,
            establishment_probability,
            generator);
        mortality_tracker = tmp.back();
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
};

}  // namespace pops

#endif  // POPS_SIMULATION_HPP
