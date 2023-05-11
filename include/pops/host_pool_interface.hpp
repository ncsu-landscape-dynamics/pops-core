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

#ifndef POPS_HOST_POOL_INTERFACE_HPP
#define POPS_HOST_POOL_INTERFACE_HPP

#include "environment_interface.hpp"

namespace pops {

template<
    typename IntegerRaster,
    typename FloatRaster,
    typename RasterIndex,
    typename Generator>
class HostPoolInterface
{
public:
    // virtual ~HostPoolInterface() = 0;
    virtual int
    disperser_to(RasterIndex row, RasterIndex col, Generator& generator) = 0;
    virtual void add_disperser_at(RasterIndex row, RasterIndex col) = 0;
    virtual double establishment_probability_at(
        RasterIndex row, RasterIndex col, IntegerRaster& susceptible) = 0;
    virtual int pest_from(RasterIndex i, RasterIndex j, int count) = 0;
    virtual int pests_to(RasterIndex row, RasterIndex col, int count) = 0;
    virtual int move_hosts_from_to(
        RasterIndex row_from,
        RasterIndex col_from,
        RasterIndex row_to,
        RasterIndex col_to,
        int count,
        Generator& generator) = 0;
    virtual void remove_infected_at(
        RasterIndex i, RasterIndex j, int count, Generator& generator) = 0;
    virtual void remove_exposed_at(
        RasterIndex i, RasterIndex j, int count, Generator& generator) = 0;

    // Brings exposed dependency to more items, needs to wait for more complete host.
    /*
    template<typename Generator>
    virtual void remove_infection_at(
        RasterIndex i,
        RasterIndex j,
        double percentage,
        IntegerRaster& infected,
        IntegerRaster& susceptible,
        std::vector<IntegerRaster>& exposed,
        IntegerRaster& total_exposed,
        Generator& generator) = 0;
        */

    // TODO: Rate and time lag will eventually go to constructor (host properties).
    virtual void apply_mortality_at(
        RasterIndex i,
        RasterIndex j,
        double mortality_rate,
        int mortality_time_lag) = 0;
    virtual int infected_at(RasterIndex i, RasterIndex j) const = 0;
    virtual int susceptible_at(RasterIndex i, RasterIndex j) const = 0;
    virtual int exposed_at(RasterIndex i, RasterIndex j) const = 0;
    virtual int total_hosts_at(RasterIndex i, RasterIndex j) const = 0;
    virtual void step_forward_mortality() = 0;
};

}  // namespace pops

#endif  // POPS_HOST_POOL_INTERFACE_HPP
