/*
 * PoPS model - Pest-host-use table for hosts and pest
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

#ifndef POPS_PEST_HOST_USE_TABLE_HPP
#define POPS_PEST_HOST_USE_TABLE_HPP

#include <vector>
#include <stdexcept>
#include <string>

namespace pops {

template<typename HostPool>
class PestHostUseTable
{
public:
    using Environment = typename HostPool::Environment;

    PestHostUseTable(const Environment& environment) : environment_(environment) {}

    void add_host_info(double susceptibility)
    {
        susceptibilities_.push_back(susceptibility);
    }

    double susceptibility(const HostPool* host) const
    {
        auto host_index = environment_.host_index(host);
        return susceptibilities_.at(host_index);
    }

private:
    std::vector<double> susceptibilities_;
    const Environment& environment_;
};

}  // namespace pops

#endif  // POPS_PEST_HOST_USE_TABLE_HPP
