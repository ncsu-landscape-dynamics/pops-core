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

    void step_forward(unsigned step)
    {
        for (auto& item : host_pools_) {
            item->step_forward(step);
        }
    }

    const std::vector<std::vector<int>>& suitable_cells() const
    {
        // TODO: if host_pools_ is empty
        return host_pools_[0]->suitable_cells();
    }

    void remove_all_infected_at(RasterIndex row, RasterIndex col, Generator& generator)
    {
        for (auto& item : host_pools_) {
            item->remove_all_infected_at(row, col, generator);
        }
    }

private:
    std::vector<HostPool*> host_pools_;  // non-owning
};

}  // namespace pops

#endif  // POPS_MULTI_HOST_POOL_HPP
