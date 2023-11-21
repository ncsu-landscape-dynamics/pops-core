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
#include "pest_host_use_table.hpp"
#include "utils.hpp"

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

    void set_pest_host_use_table(const PestHostUseTable<HostPool>& table)
    {
        for (auto& host_pool : host_pools_) {
            host_pool->set_pest_host_use_table(table);
        }
    }

    void set_competency_table(const CompetencyTable<HostPool, RasterIndex>& table)
    {
        for (auto& host_pool : host_pools_) {
            host_pool->set_competency_table(table);
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

    int total_hosts_at(RasterIndex row, RasterIndex col) const
    {
        int sum = 0;
        for (auto& host_pool : host_pools_) {
            sum += host_pool->total_hosts_at(row, col);
        }
        return sum;
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

    /**
     * @brief Move pests from a cell (multi-host)
     *
     * This is a multi-host version of single host pests_from method.
     * It randomly distributes the number to un-infect between multiple
     * hosts.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param count Number of pests requested to move from the cell
     * @param generator Random number generator
     *
     * @return Number of pests actually moved from the cell
     *
     * @note For consitency with the previous implementation, this does not modify
     * mortality cohorts nor touches the exposed cohorts.
     */
    int pests_from(RasterIndex row, RasterIndex col, int count, Generator& generator)
    {
        std::vector<int> infected;
        int index = 0;
        for (auto& host_pool : host_pools_) {
            infected.insert(infected.end(), host_pool->infected_at(row, col), index);
            index++;
        }

        index = 0;
        int collect_count = 0;
        std::vector<int> draw = draw_n_from_v(infected, count, generator);
        for (auto& host_pool : host_pools_) {
            count = std::count(draw.begin(), draw.end(), index);
            collect_count += host_pool->pests_from(row, col, count, generator);
            index++;
        }

        return collect_count;
    }
    /**
     * @brief Move pests to a cell (multi-host)
     *
     * This is a multi-host version of single host pests_to method.
     * It randomly distributes the number to infect between multiple
     * hosts.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param count Number of pests requested to move to the cell
     * @param generator Random number generator
     *
     * @return Number of accepted pests
     *
     * @note For consistency with the previous implementation, this does not make hosts
     * exposed in the SEI model. Instead, the hosts are infected right away. This may
     * become a feature in the future.
     *
     * @note For consistency with the previous implementation, this does not modify the
     * mortality cohorts. This wil need to be fixed in the future.
     *
     * @note This may be merged with add_disperser_at() in the future.
     */
    int pests_to(RasterIndex row, RasterIndex col, int count, Generator& generator)
    {
        std::vector<int> susceptible;
        int index = 0;
        for (auto& host_pool : host_pools_) {
            susceptible.insert(
                susceptible.end(), host_pool->susceptible_at(row, col), index);
            index++;
        }
        index = 0;
        int collect_count = 0;
        std::vector<int> draw = draw_n_from_v(susceptible, count, generator);
        for (auto& host_pool : host_pools_) {
            count = std::count(draw.begin(), draw.end(), index);
            collect_count += host_pool->pests_to(row, col, count, generator);
            index++;
        }

        return collect_count;
    }

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

    std::vector<HostPool*>& host_pools()
    {
        return host_pools_;
    }

private:
    std::vector<HostPool*> host_pools_;  // non-owning
};

}  // namespace pops

#endif  // POPS_MULTI_HOST_POOL_HPP
