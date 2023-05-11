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
#include "host_pool_interface.hpp"

namespace pops {

/** The type of a epidemiological model (SI or SEI)
 */
enum class ModelType
{
    SusceptibleInfected,  ///< SI (susceptible - infected)
    SusceptibleExposedInfected  ///< SEI (susceptible - exposed - infected)
};

/*! Get a corresponding enum value for a string which is a model type name.
 *
 * Throws an std::invalid_argument exception if the value was not
 * found or is not supported (which is the same thing).
 */
inline ModelType model_type_from_string(const std::string& text)
{
    if (text == "SI" || text == "SusceptibleInfected" || text == "susceptible-infected"
        || text == "susceptible_infected")
        return ModelType::SusceptibleInfected;
    else if (
        text == "SEI" || text == "SusceptibleExposedInfected"
        || text == "susceptible-exposed-infected"
        || text == "susceptible_exposed_infected")
        return ModelType::SusceptibleExposedInfected;
    else
        throw std::invalid_argument(
            "model_type_from_string: Invalid"
            " value '"
            + text + "' provided");
}

/*! Overload which allows to pass C-style string which is nullptr (NULL)
 *
 * @see model_type_from_string(const std::string& text)
 */
inline ModelType model_type_from_string(const char* text)
{
    // call the string version
    return model_type_from_string(text ? std::string(text) : std::string());
}

template<
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename Generator>
class HostPool
    : public HostPoolInterface<IntegerRaster, FloatRaster, RasterIndex, Generator>
{
public:
    using Environment1 =
        EnvironmentInterface<IntegerRaster, FloatRaster, RasterIndex, Generator>;

    HostPool(
        ModelType model_type,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& infected,
        IntegerRaster& total_exposed,
        IntegerRaster& resistant,
        std::vector<IntegerRaster>& mortality_tracker_vector,
        IntegerRaster& died,
        IntegerRaster& total_hosts,
        const Environment1& environment,
        bool establishment_stochasticity,
        double establishment_probability,
        std::vector<std::vector<int>>& suitable_cells)
        : susceptible_(susceptible),
          infected_(infected),
          exposed_(exposed),
          total_exposed_(total_exposed),
          resistant_(resistant),
          mortality_tracker_vector_(mortality_tracker_vector),
          died_(died),
          total_hosts_(total_hosts),
          environment_(environment),
          model_type_(model_type),
          establishment_stochasticity_(establishment_stochasticity),
          deterministic_establishment_probability_(establishment_probability),
          suitable_cells_(suitable_cells)
    {}

    /**
     * @brief Move (add) disperser to a cell in the host pool
     *
     * Processes event when a disperser lands in a cell potentially establishing on a
     * host. The disperser may or may not establish a based on host availability,
     * weather, establishment probability, and stochasticity.
     *
     * Any new dispersers targeting host in the host pool should be added using this
     * function.
     *
     * @param row Row number of the target cell
     * @param col Column number of the target cell
     * @param susceptible Raster of susceptible hosts
     * @param exposed_or_infected Raster of exposed or infected hosts
     * @param mortality_tracker Raster tracking hosts for mortality
     * @param total_populations Raster of all individuals (hosts and non-hosts)
     * @param total_exposed Raster tracking all exposed hosts
     * @param weather Whether to use weather
     * @param weather_coefficient Raster with weather coefficients per cell
     * @param establishment_probability Fixed probability disperser establishment
     *
     * @return true if disperser has established in the cell, false otherwise
     *
     * @throw std::runtime_error if model type is unsupported (i.e., not SI or SEI)
     */
    int disperser_to(
        RasterIndex row,
        RasterIndex col,
        IntegerRaster& mortality_tracker,
        const IntegerRaster& total_populations,
        Generator& generator)
    {
        std::uniform_real_distribution<double> distribution_uniform(0.0, 1.0);
        if (susceptible_(row, col) > 0) {
            double probability_of_establishment =
                establishment_probability_at(row, col, susceptible_, total_populations);
            double establishment_tester = 1 - deterministic_establishment_probability_;
            if (establishment_stochasticity_)
                establishment_tester = distribution_uniform(generator);
            if (establishment_tester < probability_of_establishment) {
                add_disperser_at(row, col);
                susceptible_(row, col) -= 1;
                if (model_type_ == ModelType::SusceptibleInfected) {
                    mortality_tracker(row, col) += 1;
                }
                else if (model_type_ == ModelType::SusceptibleExposedInfected) {
                    total_exposed_(row, col) += 1;
                }
                else {
                    throw std::runtime_error(
                        "Unknown ModelType value in "
                        "Simulation::disperse()");
                }
                return true;
            }
        }
        return false;
    }

    void add_disperser_at(RasterIndex i, RasterIndex j)
    {
        if (model_type_ == ModelType::SusceptibleInfected) {
            infected_(i, j) += 1;
        }
        else if (model_type_ == ModelType::SusceptibleExposedInfected) {
            exposed_.back()(i, j) += 1;
        }
        else {
            throw std::runtime_error(
                "Unknown ModelType value in HostPool::add_disperser_at()");
        }
    }

    double establishment_probability_at(
        RasterIndex row,
        RasterIndex col,
        IntegerRaster& susceptible,
        const IntegerRaster& total_populations)
    {
        double probability_of_establishment =
            (double)(susceptible(row, col)) / total_populations(row, col);
        return environment_.influence_probability_of_establishment_at(
            row, col, probability_of_establishment);
    }

    int pest_from(RasterIndex i, RasterIndex j, int count)
    {
        susceptible_(i, j) += count;
        infected_(i, j) -= count;
        return count;
    }

    int pests_to(RasterIndex row, RasterIndex col, int count)
    {
        // The target cell can accept all.
        if (susceptible_(row, col) >= count) {
            susceptible_(row, col) -= count;
            infected_(row, col) += count;
        }
        // More pests than the target cell can accept.
        // This can happen if there is simply not enough S hosts to accommodate all
        // the pests moving from the source or if multiple sources end up in the
        // same target cell and there is not enough S hosts to accommodate all of
        // them. The pests just disappear in both cases.
        else {
            count = susceptible_(row, col);
            susceptible_(row, col) -= count;
            infected_(row, col) += count;
        }
        return count;
    }

    int move_hosts_from_to(
        RasterIndex row_from,
        RasterIndex col_from,
        RasterIndex row_to,
        RasterIndex col_to,
        int count,
        Generator& generator)
    {
        int total_hosts_moved{count};
        if (count > total_hosts_(row_from, col_from)) {
            total_hosts_moved = total_hosts_(row_from, col_from);
        }
        int total_infecteds = infected_(row_from, col_from);
        int suscepts = susceptible_(row_from, col_from);
        int expose = total_exposed_(row_from, col_from);
        int resist = resistant_(row_from, col_from);
        // set up vector of numeric categories (infected = 1, susceptible = 2,
        // exposed = 3) for drawing # moved in each category
        std::vector<int> categories(total_infecteds, 1);
        categories.insert(categories.end(), suscepts, 2);
        categories.insert(categories.end(), expose, 3);
        categories.insert(categories.end(), resist, 4);

        std::vector<int> draw = draw_n_from_v(categories, total_hosts_moved, generator);
        int infected_moved = std::count(draw.begin(), draw.end(), 1);
        int susceptible_moved = std::count(draw.begin(), draw.end(), 2);
        int exposed_moved = std::count(draw.begin(), draw.end(), 3);
        int resistant_moved = std::count(draw.begin(), draw.end(), 4);

        if (exposed_moved > 0) {
            std::vector<int> exposed_draw = draw_n_from_cohorts(
                exposed_, exposed_moved, row_from, col_from, generator);
            int index = 0;
            for (auto& raster : exposed_) {
                raster(row_from, col_from) -= exposed_draw[index];
                raster(row_to, col_to) += exposed_draw[index];
                index += 1;
            }
        }
        if (infected_moved > 0) {
            std::vector<int> mortality_draw = draw_n_from_cohorts(
                mortality_tracker_vector_,
                infected_moved,
                row_from,
                col_from,
                generator);
            int index = 0;
            for (auto& raster : mortality_tracker_vector_) {
                raster(row_from, col_from) -= mortality_draw[index];
                raster(row_to, col_to) += mortality_draw[index];
                index += 1;
            }
        }
        // check that the location with host movement is in suitable cells.
        // Since suitable-cells comes from the total hosts originally. The
        // the first check is for total_hosts
        if (total_hosts_(row_to, col_to) == 0) {
            for (auto indices : suitable_cells_) {
                int i = indices[0];
                int j = indices[1];
                if ((i == row_to) && (j == col_to)) {
                    std::vector<int> added_index = {row_to, col_to};
                    suitable_cells_.push_back(added_index);
                    break;
                }
            }
        }

        infected_(row_from, col_from) -= infected_moved;
        susceptible_(row_from, col_from) -= susceptible_moved;
        total_hosts_(row_from, col_from) -= total_hosts_moved;
        total_exposed_(row_from, col_from) -= exposed_moved;
        resistant_(row_from, col_from) -= resistant_moved;
        infected_(row_to, col_to) += infected_moved;
        susceptible_(row_to, col_to) += susceptible_moved;
        total_hosts_(row_to, col_to) += total_hosts_moved;
        total_exposed_(row_to, col_to) += exposed_moved;
        resistant_(row_to, col_to) += resistant_moved;
        return total_hosts_moved;
    }

    void
    remove_infected_at(RasterIndex i, RasterIndex j, int count, Generator& generator)
    {
        // remove percentage of infestation/infection in the infected class
        infected_(i, j) -= count;
        // remove the removed infected from mortality cohorts
        if (count > 0) {
            std::vector<int> mortality_draw =
                draw_n_from_cohorts(mortality_tracker_vector_, count, i, j, generator);
            int index = 0;
            for (auto& raster : mortality_tracker_vector_) {
                raster(i, j) -= mortality_draw[index];
                index += 1;
            }
        }
        // move infested/infected host back to susceptible pool
        susceptible_(i, j) += count;
    }

    void
    remove_exposed_at(RasterIndex i, RasterIndex j, int count, Generator& generator)
    {
        // remove the same percentage for total exposed and remove randomly from
        // each cohort
        total_exposed_(i, j) -= count;
        if (count > 0) {
            std::vector<int> exposed_draw =
                draw_n_from_cohorts(exposed_, count, i, j, generator);
            int index = 0;
            for (auto& raster : exposed_) {
                raster(i, j) -= exposed_draw[index];
                index += 1;
            }
        }
        // move infested/infected host back to susceptible pool
        susceptible_(i, j) += count;
    }

    // Brings exposed dependency to more items, needs to wait for more complete host.
    /*
    template<typename Generator>
    void remove_infection_at(
        RasterIndex i,
        RasterIndex j,
        double percentage,
        IntegerRaster& infected,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& total_exposed,
        Generator& generator)
    {
        this->remove_infected_at(i, j, 1, infected, susceptible, generator);
        this->remove_exposed_at(i, j, 1, susceptible, exposed, generator);
    }
    */

    // TODO: Rate and time lag will eventually go to constructor (host properties).
    void apply_mortality_at(
        RasterIndex i, RasterIndex j, double mortality_rate, int mortality_time_lag)
    {
        int max_index = mortality_tracker_vector_.size() - mortality_time_lag - 1;
        for (int index = 0; index <= max_index; index++) {
            int mortality_in_index = 0;
            if (mortality_tracker_vector_[index](i, j) > 0) {
                // used to ensure that all infected hosts in the last year of
                // tracking mortality
                if (index == 0) {
                    mortality_in_index = mortality_tracker_vector_[index](i, j);
                }
                else {
                    mortality_in_index =
                        mortality_rate * mortality_tracker_vector_[index](i, j);
                }
                mortality_tracker_vector_[index](i, j) -= mortality_in_index;
                died_(i, j) += mortality_in_index;
                if (infected_(i, j) > 0) {
                    infected_(i, j) -= mortality_in_index;
                }
                if (total_hosts_(i, j) > 0) {
                    total_hosts_(i, j) -= mortality_in_index;
                }
            }
        }
    }

    int infected_at(RasterIndex i, RasterIndex j)
    {
        return infected_(i, j);
    }

    int susceptible_at(RasterIndex i, RasterIndex j)
    {
        return susceptible_(i, j);
    }

    int exposed_at(RasterIndex i, RasterIndex j)
    {
        // Future code could remove total exposed and compute that on the fly.
        //        int sum = 0;
        //        for (const auto& raster : exposed_)
        //            sum += raster(i, j);
        //        return sum;
        return total_exposed_(i, j);
    }

    int total_hosts_at(RasterIndex i, RasterIndex j)
    {
        // computed instead of using a raster
        return susceptible_at(i, j) /*+ exposed_at(i, j)*/ + infected_at(i, j);
    }

    void step_forward_mortality()
    {
        rotate_left_by_one(mortality_tracker_vector_);
    }

private:
    IntegerRaster& susceptible_;
    IntegerRaster& infected_;

    std::vector<IntegerRaster>& exposed_;
    IntegerRaster& total_exposed_;

    IntegerRaster& resistant_;

    std::vector<IntegerRaster>& mortality_tracker_vector_;
    IntegerRaster& died_;

    IntegerRaster& total_hosts_;
    const Environment1& environment_;

    ModelType model_type_;

    bool establishment_stochasticity_{true};
    double deterministic_establishment_probability_{0};

    std::vector<std::vector<int>>& suitable_cells_;
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
                hosts.remove_infected_at(i, j, removed_infected, generator);
                // remove the same percentage for total exposed and remove randomly from
                // each cohort
                auto exposed = hosts.exposed_at(i, j);
                int total_removed_exposed =
                    exposed - std::lround(exposed * survival_rate_(i, j));
                hosts.remove_exposed_at(i, j, total_removed_exposed, generator);
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
                hosts.remove_infected_at(i, j, count, generator);
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
                row_from, col_from, row_to, col_to, moved[4], generator);
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
    typename Generator = std::default_random_engine>
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
    Generator generator_;
    /// Non-owning pointer to environment for weather
    const Environment<IntegerRaster, FloatRaster, RasterIndex>* environment_{nullptr};
    /**
     * Optional soil pool
     */
    std::shared_ptr<SoilPool<IntegerRaster, FloatRaster, RasterIndex>> soil_pool_{
        nullptr};
    /**
     * Percentage (0-1 ratio) of disperers to be send to soil
     */
    double to_soil_percentage_{0};

public:
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
        unsigned random_seed,
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
    {
        generator_.seed(random_seed);
    }

    Simulation() = delete;

    /**
     * @brief Set environment used for weather to provided environment
     * @param environment Pointer to an existing environment
     *
     * The simulation object does not take ownership of the environment.
     */
    void set_environment(
        const Environment<IntegerRaster, FloatRaster, RasterIndex>* environment)
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
    const Environment<IntegerRaster, FloatRaster, RasterIndex>*
    environment(bool allow_empty = false)
    {
        static Environment<IntegerRaster, FloatRaster, RasterIndex> empty;
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
    void remove(
        IntegerRaster& infected,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& total_exposed,
        std::vector<IntegerRaster>& mortality_tracker_vector,
        const FloatRaster& temperature,
        double lethal_temperature,
        std::vector<std::vector<int>>& suitable_cells)
    {
        IntegerRaster empty;
        StandardHostPool hosts(
            model_type_,
            infected,
            exposed,
            susceptible,
            total_exposed,
            empty,
            mortality_tracker_vector,
            empty,
            empty,
            *environment(true),
            false,
            0,
            suitable_cells);
        RemoveByTemperature<StandardHostPool, IntegerRaster, FloatRaster> remove;
        remove.action(
            hosts, temperature, lethal_temperature, suitable_cells, generator_);
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
    void remove_percentage(
        IntegerRaster& infected,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& mortality_tracker_vector,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& total_exposed,
        const FloatRaster& survival_rate,
        std::vector<std::vector<int>>& suitable_cells)
    {
        IntegerRaster empty;
        StandardHostPool hosts(
            model_type_,
            susceptible,
            exposed,
            infected,
            total_exposed,
            empty,
            mortality_tracker_vector,
            empty,
            empty,
            *environment(true),
            false,
            0,
            suitable_cells);
        SurvivalRateAction<StandardHostPool, IntegerRaster, FloatRaster> survival(
            survival_rate);
        survival.action(hosts, suitable_cells, generator_);
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
            infected,
            empty,
            empty,
            mortality_tracker_vector,
            died,
            total_hosts,
            *environment(true),
            false,
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
        std::vector<std::vector<int>>& suitable_cells)
    {
        HostMovement<StandardHostPool, IntegerRaster, FloatRaster, RasterIndex>
            host_movement{};
        IntegerRaster empty;
        StandardHostPool hosts{
            model_type_,
            susceptible,
            exposed,
            infected,
            total_exposed,
            resistant,
            mortality_tracker_vector,
            empty,
            total_hosts,
            *environment(true),
            false,
            0,
            suitable_cells};
        return host_movement.movement(
            hosts, step, last_index, movements, movement_schedule, generator_);
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
        const std::vector<std::vector<int>>& suitable_cells)
    {
        double lambda = reproductive_rate;
        for (auto indices : suitable_cells) {
            int i = indices[0];
            int j = indices[1];
            if (infected(i, j) > 0) {
                if (weather)
                    lambda =
                        reproductive_rate * environment()->weather_coefficient_at(i, j);
                int dispersers_from_cell = 0;
                if (dispersers_stochasticity_) {
                    std::poisson_distribution<int> distribution(lambda);
                    for (int k = 0; k < infected(i, j); k++) {
                        dispersers_from_cell += distribution(generator_);
                    }
                }
                else {
                    dispersers_from_cell = lambda * infected(i, j);
                }
                if (soil_pool_) {
                    // From all the generated dispersers, some go to the soil in the
                    // same cell and don't participate in the kernel-driven dispersal.
                    auto dispersers_to_soil =
                        std::round(to_soil_percentage_ * dispersers_from_cell);
                    soil_pool_->dispersers_to(dispersers_to_soil, i, j, generator_);
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
        double establishment_probability = 0.5)
    {
        // The interaction does not happen over the member variables yet, use empty
        // variables. This requires SI/SEI to be fully resolved in host and not in
        // disperse_and_infect.
        IntegerRaster empty;
        std::vector<IntegerRaster> empty_vector;
        StandardHostPool host_pool{
            model_type_,
            susceptible,
            exposed,
            infected,
            total_exposed,
            empty,
            empty_vector,
            empty,
            empty,
            *environment(!weather),
            establishment_stochasticity_,
            establishment_probability,
            suitable_cells};

        std::uniform_real_distribution<double> distribution_uniform(0.0, 1.0);
        int row;
        int col;

        for (auto indices : suitable_cells) {
            int i = indices[0];
            int j = indices[1];
            if (dispersers(i, j) > 0) {
                for (int k = 0; k < dispersers(i, j); k++) {
                    std::tie(row, col) = dispersal_kernel(generator_, i, j);
                    if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
                        // export dispersers dispersed outside of modeled area
                        outside_dispersers.emplace_back(std::make_tuple(row, col));
                        established_dispersers(i, j) -= 1;
                        continue;
                    }
                    // Put a disperser to the host pool.
                    auto dispersed = host_pool.disperser_to(
                        row, col, mortality_tracker, total_populations, generator_);
                    if (!dispersed) {
                        established_dispersers(i, j) -= 1;
                    }
                }
            }
            if (soil_pool_) {
                // Get dispersers from the soil if there are any.
                auto num_dispersers = soil_pool_->dispersers_from(i, j, generator_);
                // Put each disperser to the host pool.
                for (int k = 0; k < num_dispersers; k++) {
                    host_pool.disperser_to(
                        i, j, mortality_tracker, total_populations, generator_);
                }
            }
        }
    }

    // For backwards compatibility for tests.
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
        std::vector<IntegerRaster> empty_vector;
        disperse(
            dispersers,
            established_dispersers,
            susceptible,
            empty_vector,
            infected,
            mortality_tracker,
            total_populations,
            total_exposed,
            outside_dispersers,
            weather,
            dispersal_kernel,
            suitable_cells,
            establishment_probability);
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
        double leaving_percentage)
    {
        UNUSED(total_hosts);  // Total hosts is computed now.
        IntegerRaster empty;
        std::vector<IntegerRaster> empty_vector;
        StandardHostPool hosts{
            model_type_,
            susceptible,
            empty_vector,
            infected,
            empty,
            empty,
            empty_vector,
            empty,
            empty,
            *environment(true),
            false,
            0,
            suitable_cells};
        MoveOverpopulatedPests<
            StandardHostPool,
            IntegerRaster,
            FloatRaster,
            RasterIndex>
            move_pest{overpopulation_percentage, leaving_percentage, rows_, cols_};
        move_pest.action(
            hosts, outside_dispersers, dispersal_kernel, suitable_cells, generator_);
    }

    /** Infect exposed hosts (E to I transition in the SEI model)
     *
     * Applicable to SEI model, no-operation otherwise, i.e., parameters
     * are left intact for other models.
     *
     * The exposed vector are the hosts exposed in the previous steps.
     * The length of the vector is the number of steps of the latency
     * period plus one. Before the first latency period is over,
     * the E to I transition won't happen because no item in the exposed
     * vector is old enough to become infected.
     *
     * The position of the items in the exposed vector determines their
     * age, i.e., for how long the hosts are exposed. The oldest item
     * is at the front and youngest at the end.
     * Before the the first latency period is over, items in the front
     * are still empty (unused) because no hosts were exposed for the
     * given time period.
     * After the first latency
     * period, this needs to be true before the function is called and
     * it is true after the function
     * finished with the difference that after the function is called,
     * the last item is empty in the sense that it does not contain any
     * hosts.
     *
     * When the E to I transition happens, hosts from the oldest item
     * in the exposed vector are moved to the infected (and mortality
     * tracker). They are removed from the exposed item and this item
     * is moved to the back of the vector.
     *
     * Like in disperse(), there is no distinction between *infected*
     * and *mortality_tracker*, but different usage is expected outside
     * of this function.
     *
     * The raster class used with the simulation class needs to support
     * `.fill()` method for this function to work.
     *
     * @param step Step in the simulation (>=0)
     * @param exposed Vector of exposed hosts
     * @param infected Infected hosts
     * @param mortality_tracker Newly infected hosts
     * @param total_exposed Total exposed in all exposed cohorts
     */
    void infect_exposed(
        unsigned step,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& infected,
        IntegerRaster& mortality_tracker,
        IntegerRaster& total_exposed)
    {
        if (model_type_ == ModelType::SusceptibleExposedInfected) {
            if (step >= latency_period_) {
                // Oldest item needs to be in the front
                auto& oldest = exposed.front();
                // Move hosts to infected raster
                infected += oldest;
                mortality_tracker += oldest;
                total_exposed += (oldest * (-1));
                // Reset the raster
                // (hosts moved from the raster)
                oldest.fill(0);
            }
            // Age the items and the used one to the back
            // elements go one position to the left
            // new oldest goes to the front
            // old oldest goes to the back
            rotate_left_by_one(exposed);
        }
        else if (model_type_ == ModelType::SusceptibleInfected) {
            // no-op
        }
        else {
            throw std::runtime_error(
                "Unknown ModelType value in Simulation::infect_exposed()");
        }
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
        IntegerRaster& mortality_tracker,
        const IntegerRaster& total_populations,
        IntegerRaster& total_exposed,
        std::vector<std::tuple<int, int>>& outside_dispersers,
        bool weather,
        DispersalKernel& dispersal_kernel,
        std::vector<std::vector<int>>& suitable_cells,
        double establishment_probability = 0.5)
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
            establishment_probability);
        if (model_type_ == ModelType::SusceptibleExposedInfected) {
            this->infect_exposed(
                step, exposed, infected, mortality_tracker, total_exposed);
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
        std::shared_ptr<SoilPool<IntegerRaster, FloatRaster, RasterIndex>> soil_pool,
        double dispersers_percentage)
    {
        this->soil_pool_ = soil_pool;
        this->to_soil_percentage_ = dispersers_percentage;
    }

    Generator& random_number_generator()
    {
        return generator_;
    }
};

}  // namespace pops

#endif  // POPS_SIMULATION_HPP
