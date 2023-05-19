/*
 * PoPS model - pest or pathogen spread simulation
 *
 * Copyright (C) 2023 by the authors.
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

#ifndef POPS_HOST_POOL_HPP
#define POPS_HOST_POOL_HPP

#include <vector>
#include <random>
#include <stdexcept>
#include <algorithm>

#include "host_pool_interface.hpp"
#include "model_type.hpp"
#include "environment_interface.hpp"

namespace pops {

template<
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename Generator>
class HostPool
    : public HostPoolInterface<IntegerRaster, FloatRaster, RasterIndex, Generator>
{
public:
    using Environment =
        EnvironmentInterface<IntegerRaster, FloatRaster, RasterIndex, Generator>;

    HostPool(
        ModelType model_type,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& exposed,
        unsigned latency_period,
        IntegerRaster& infected,
        IntegerRaster& total_exposed,
        IntegerRaster& resistant,
        std::vector<IntegerRaster>& mortality_tracker_vector,
        IntegerRaster& died,
        IntegerRaster& total_hosts,
        const Environment& environment,
        bool dispersers_stochasticity,
        double reproductive_rate,
        bool establishment_stochasticity,
        double establishment_probability,
        RasterIndex rows,
        RasterIndex cols,
        std::vector<std::vector<int>>& suitable_cells)
        : susceptible_(susceptible),
          infected_(infected),
          exposed_(exposed),
          latency_period_(latency_period),
          total_exposed_(total_exposed),
          resistant_(resistant),
          mortality_tracker_vector_(mortality_tracker_vector),
          died_(died),
          total_hosts_(total_hosts),
          environment_(environment),
          model_type_(model_type),
          dispersers_stochasticity_(dispersers_stochasticity),
          reproductive_rate_(reproductive_rate),
          establishment_stochasticity_(establishment_stochasticity),
          deterministic_establishment_probability_(establishment_probability),
          rows_(rows),
          cols_(cols),
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
    int disperser_to(RasterIndex row, RasterIndex col, Generator& generator)
    {
        std::uniform_real_distribution<double> distribution_uniform(0.0, 1.0);
        if (susceptible_(row, col) > 0) {
            double probability_of_establishment =
                establishment_probability_at(row, col, susceptible_);
            double establishment_tester = 1 - deterministic_establishment_probability_;
            if (establishment_stochasticity_)
                establishment_tester = distribution_uniform(generator);
            if (establishment_tester < probability_of_establishment) {
                add_disperser_at(row, col);
                return true;
            }
        }
        return false;
    }

    void add_disperser_at(RasterIndex row, RasterIndex col)
    {
        susceptible_(row, col) -= 1;
        if (model_type_ == ModelType::SusceptibleInfected) {
            infected_(row, col) += 1;
            mortality_tracker_vector_.back()(row, col) += 1;
        }
        else if (model_type_ == ModelType::SusceptibleExposedInfected) {
            exposed_.back()(row, col) += 1;
            total_exposed_(row, col) += 1;
        }
        else {
            throw std::runtime_error(
                "Unknown ModelType value in HostPool::add_disperser_at()");
        }
    }

    int dispersers_from(RasterIndex row, RasterIndex col, Generator& generator) const
    {
        if (infected_at(row, col) <= 0)
            return 0;
        double lambda =
            environment_.influence_reproductive_rate_at(row, col, reproductive_rate_);
        int dispersers_from_cell = 0;
        if (dispersers_stochasticity_) {
            std::poisson_distribution<int> distribution(lambda);
            for (int k = 0; k < infected_at(row, col); k++) {
                dispersers_from_cell += distribution(generator);
            }
        }
        else {
            dispersers_from_cell = lambda * infected_at(row, col);
        }
        return dispersers_from_cell;
    }

    double establishment_probability_at(
        RasterIndex row, RasterIndex col, IntegerRaster& susceptible)
    {
        double probability_of_establishment =
            (double)(susceptible(row, col))
            / environment_.total_population_at(row, col);
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

    // Notably, this does not remove resistant.
    void completely_remove_susceptible_at(RasterIndex row, RasterIndex col, int count)
    {
        if (count <= 0)
            return;
        susceptible_(row, col) -= count;
        reset_total_host(row, col);
    }

    // special overload for double, so that handling of floating point values is managed
    // here
    void
    completely_remove_susceptible_at(RasterIndex row, RasterIndex col, double count)
    {
        if (count <= 0)
            return;
        susceptible_(row, col) = susceptible_(row, col) - count;
        reset_total_host(row, col);
    }

    void completely_remove_exposed_at(
        RasterIndex row, RasterIndex col, std::vector<double> counts)
    {
        int total = 0;
        // no simple zip in C++, falling back to indices
        // TODO: check sizes in any case
        for (size_t i = 0; i < counts.size(); ++i) {
            exposed_[i](row, col) -= counts[i];
            total += counts[i];
        }
        reset_total_host(row, col);
    }

    // Exposed may need an additonal method which uses RNG to distrubute
    // individuals over the cohorts.

    void completely_remove_infected_at(RasterIndex row, RasterIndex col, int count)
    {
        // Possibly reuse in the I->S removal.
        if (count <= 0)
            return;
        infected_(row, col) -= count;
        // TODO: The distribution among moratlity cohorts needs to be done
        // either by passing a generator and drawing from cohorts or by
        // the caller passing a vector of how they should be distributed like
        // for exposed.
        //        std::default_random_engine generator;
        //        std::vector<int> mortality_draw =
        //            draw_n_from_cohorts(mortality_tracker_vector_, count, row, col,
        //            generator);
        //        int index = 0;
        //        for (auto& raster : mortality_tracker_vector_) {
        //            raster(row, col) -= mortality_draw[index];
        //            index += 1;
        //        }
        reset_total_host(row, col);
    }

    void completely_remove_infected_at(RasterIndex row, RasterIndex col, double count)
    {
        // Possibly reuse in the I->S removal.
        if (count <= 0)
            return;
        infected_(row, col) -= count;
        reset_total_host(row, col);
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

    void make_resistant_at(
        RasterIndex row,
        RasterIndex col,
        int susceptible,
        const std::vector<int>& exposed,
        int infected)
    {
        int total_resistant = 0;
        // TODO: check negative numbers in these cases?
        susceptible_(row, col) -= susceptible;
        total_resistant += susceptible;
        // no simple zip in C++, falling back to indices
        // TODO: check sizes in any case
        for (size_t i = 0; i < exposed.size(); ++i) {
            exposed_[i](row, col) -= exposed[i];
            total_resistant += exposed[i];
        }
        infected_(row, col) -= infected;
        // TODO: mortality cohorts need to be reduced when infected is reduced
        //        for (auto& raster : mortality_tracker_vector_) {
        //            raster(row, col) -= infected /
        //            mortality_tracker_vector_.size();
        //        }
        total_resistant += infected;
        resistant_(row, col) += total_resistant;
    }

    void remove_resistance_at(RasterIndex row, RasterIndex col)
    {
        susceptible_(row, col) += resistant_(row, col);
        resistant_(row, col) = 0;
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

    int infected_at(RasterIndex i, RasterIndex j) const
    {
        return infected_(i, j);
    }

    int susceptible_at(RasterIndex i, RasterIndex j) const
    {
        return susceptible_(i, j);
    }

    int exposed_at(RasterIndex i, RasterIndex j) const
    {
        // Future code could remove total exposed and compute that on the fly.
        return total_exposed_(i, j);
    }

    int computed_exposed_at(RasterIndex i, RasterIndex j) const
    {
        int sum = 0;
        for (const auto& raster : exposed_)
            sum += raster(i, j);
        return sum;
    }

    std::vector<int> exposed_by_group_at(RasterIndex row, RasterIndex col) const
    {
        std::vector<int> all;
        all.reserve(exposed_.size());
        for (const auto& raster : exposed_)
            all.push_back(raster(row, col));
        return all;
    }

    int resistant_at(RasterIndex row, RasterIndex col) const
    {
        return resistant_(row, col);
    }

    int total_hosts_at(RasterIndex i, RasterIndex j) const
    {
        // computed instead of using a raster
        return susceptible_at(i, j) /*+ exposed_at(i, j)*/ + infected_at(i, j);
    }

    void step_forward_mortality()
    {
        rotate_left_by_one(mortality_tracker_vector_);
    }

    bool is_outside(RasterIndex row, RasterIndex col)
    {
        return row < 0 || row >= rows_ || col < 0 || col >= cols_;
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
    void step_forward(unsigned step)
    {
        if (model_type_ == ModelType::SusceptibleExposedInfected) {
            if (step >= latency_period_) {
                // Oldest item needs to be in the front
                auto& oldest = exposed_.front();
                // Move hosts to infected raster
                infected_ += oldest;
                mortality_tracker_vector_.back() += oldest;
                total_exposed_ += (oldest * (-1));
                // Reset the raster
                // (hosts moved from the raster)
                oldest.fill(0);
            }
            // Age the items and the used one to the back
            // elements go one position to the left
            // new oldest goes to the front
            // old oldest goes to the back
            rotate_left_by_one(exposed_);
        }
        else if (model_type_ == ModelType::SusceptibleInfected) {
            // no-op
        }
        else {
            throw std::runtime_error(
                "Unknown ModelType value in Simulation::infect_exposed()");
        }
    }

private:
    // This would ideally be updated piece by piece or just computed on the fly always.
    void reset_total_host(RasterIndex row, RasterIndex col)
    {
        total_hosts_(row, col) = susceptible_(row, col) + computed_exposed_at(row, col)
                                 + infected_(row, col) + resistant_(row, col);
    }

    IntegerRaster& susceptible_;
    IntegerRaster& infected_;

    std::vector<IntegerRaster>& exposed_;
    unsigned latency_period_{0};
    IntegerRaster& total_exposed_;

    IntegerRaster& resistant_;

    std::vector<IntegerRaster>& mortality_tracker_vector_;
    IntegerRaster& died_;

    IntegerRaster& total_hosts_;
    const Environment& environment_;

    ModelType model_type_;

    bool dispersers_stochasticity_{false};
    double reproductive_rate_{0};
    bool establishment_stochasticity_{true};
    double deterministic_establishment_probability_{0};

    RasterIndex rows_{0};
    RasterIndex cols_{0};

    std::vector<std::vector<int>>& suitable_cells_;
};

}  // namespace pops

#endif  // POPS_HOST_POOL_HPP
