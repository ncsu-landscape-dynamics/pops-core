/*
 * PoPS model - random uniform dispersal kernel
 *
 * Copyright (C) 2015-2020 by the authors.
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

#ifndef POPS_NETWORK_KERNEL_HPP
#define POPS_NETWORK_KERNEL_HPP

#include "kernel_types.hpp"

#include <random>

namespace pops {

class Network
{
protected:
public:
    template<typename Generator>
    std::tuple<int, int>
    time_to_row_col(int row, int col, double time, Generator& generator)
    {
        auto node_id = get_node(row, col);
        next_node(node_id, generator);
    }
};

/*! Dispersal kernel for random uniform dispersal over the whole
 * landscape
 *
 * This class is a good example of how to write a kernel and
 * it is useful for testing due to its simplicity. It tends to generate
 * a lot of spread because it quickly spreads over the landscape.
 * However, it may work as a good starting point for cases where no
 * theory about the spread is available.
 */
class NetworkDispersalKernel
{
protected:
    Network network_;
    std::uniform_real_distribution<> time_distribution;

public:
    NetworkDispersalKernel(Network& network)
        : row_max_(row_max), col_max_(col_max), time_distribution(0, row_max)
    {}

    /*! \copybrief RadialDispersalKernel::operator()()
     *
     */
    template<typename Generator>
    std::tuple<int, int> operator()(Generator& generator, int row, int col)
    {
        double time = time_distribution(generator);
        std::tie(row, col) = network_.time_to_row_col(row, col, time, generator);

        return std::make_tuple(row, col);
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        return type == DispersalKernelType::Network;
    }
};

}  // namespace pops

#endif  // POPS_NETWORK_KERNEL_HPP
