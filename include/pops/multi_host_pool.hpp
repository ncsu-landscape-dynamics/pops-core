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

#ifndef POPS_MULTI_HOST_POOL_HPP
#define POPS_MULTI_HOST_POOL_HPP

#include <vector>
#include <algorithm>
#include <random>

#include "competency_table.hpp"

namespace pops {

template<
    typename HostPool,
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename GeneratorProvider>
class MultiHostPool
{
public:
    /**
     * Standard random number generator to be passed directly to the methods.
     */
    using Generator = typename GeneratorProvider::Generator;

    MultiHostPool(const std::vector<HostPool*>& host_pools) : host_pools_(host_pools) {}

    void
    set_competency_table(const CompetencyTable<HostPool, RasterIndex>& competency_table)
    {
        for (auto& host_pool : host_pools_) {
            host_pool->set_competency_table(competency_table);
        }
    }

    const std::vector<std::vector<int>>& suitable_cells() const
    {
        // TODO: if host_pools_ is empty
        return host_pools_[0]->suitable_cells();
    }

    bool is_outside(RasterIndex row, RasterIndex col)
    {
        // TODO: if host_pools_ is empty
        return host_pools_[0]->is_outside(row, col);
    }

    void step_forward(unsigned step)
    {
        for (auto& host_pool : host_pools_) {
            host_pool->step_forward(step);
        }
    }

    void remove_all_infected_at(RasterIndex row, RasterIndex col, Generator& generator)
    {
        for (auto& host_pool : host_pools_) {
            host_pool->remove_all_infected_at(row, col, generator);
        }
    }

    void remove_infection_by_ratio_at(
        RasterIndex row, RasterIndex col, double ratio, Generator& generator)
    {
        for (auto& host_pool : host_pools_) {
            host_pool->remove_infection_by_ratio_at(row, col, ratio, generator);
        }
    }

    void apply_mortality_at(
        RasterIndex row, RasterIndex col, double mortality_rate, int mortality_time_lag)
    {
        for (auto& host_pool : host_pools_) {
            host_pool->apply_mortality_at(row, col, mortality_rate, mortality_time_lag);
        }
    }

    void step_forward_mortality()
    {
        for (auto& host_pool : host_pools_) {
            host_pool->step_forward_mortality();
        }
    }

    /**
     * @brief Get number of infected hosts at a given cell
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     *
     * @return Number of infected hosts
     */
    int infected_at(RasterIndex row, RasterIndex col) const
    {
        int infected = 0;
        for (auto& host_pool : host_pools_) {
            infected += host_pool->infected_at(row, col);
        }
        return infected;
    }

    bool do_establishment_test(double value)
    {
        UNUSED(value);
        return true;
    }

    HostPool* pick_host_by_probability(
        std::vector<HostPool*>& hosts,
        const std::vector<double>& probabilities,
        Generator& generator)
    {
        std::discrete_distribution<int> distribution{
            probabilities.begin(), probabilities.end()};
        return hosts.at(distribution(generator));
    }

    int dispersers_from(RasterIndex row, RasterIndex col, Generator& generator) const
    {
        int sum{0};
        for (auto& host_pool : host_pools_) {
            sum += host_pool->dispersers_from(row, col, generator);
        }
        return sum;
    }

    template<typename Generator>
    int disperser_to(RasterIndex row, RasterIndex col, Generator& generator)
    {
        std::vector<HostPool*> hosts;
        std::vector<double> probabilities;
        double total_s_score = 0;
        for (auto& host_pool : host_pools_) {
            double s_for_item =
                host_pool->establishment_probability_at(row, col)  // accounts for W, N
                * host_pool->susceptibility();
            hosts.push_back(host_pool);
            probabilities.push_back(s_for_item);
            total_s_score += s_for_item;  // we should make sure this is <=1
        }
        std::string pest_or_pathogen = "pathogen";
        if (pest_or_pathogen == "pest") {
            if (do_establishment_test(total_s_score)) {  // this is now only in host
                auto host = pick_host_by_probability(hosts, probabilities, generator);
                host->add_disperser_at(row, col);  // simply increases the counts
                return 1;
            }
        }
        else if (pest_or_pathogen == "pathogen") {
            auto host = pick_host_by_probability(hosts, probabilities, generator);
            return host->disperser_to(row, col, generator);  // with establishment test
        }
        throw std::invalid_argument(
            "Unknown value for pest_or_pathogen: " + pest_or_pathogen);
    }

private:
    std::vector<HostPool*> host_pools_;  // non-owning
};

}  // namespace pops

#endif  // POPS_MULTI_HOST_POOL_HPP
