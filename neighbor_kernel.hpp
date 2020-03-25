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

#ifndef POPS_NEIGHBOR_KERNEL_HPP
#define POPS_NEIGHBOR_KERNEL_HPP

#include "kernel_types.hpp"
#include "radial_kernel.hpp"

#include <random>

namespace pops {

/*! Dispersal kernel for random uniform dispersal over the whole
 * landscape
 *
 * This class is a good example of how to write a kernel and
 * it is useful for testing due to its simplicity. It tends to generate
 * a lot of spread because it quickly spreads over the landscape.
 * However, it may work as a good starting point for cases where no
 * theory about the spread is available.
 */
class NeighborDispersalKernel
{
protected:
    Direction direction_;
public:
    NeighborDispersalKernel(Direction dispersal_direction)
        :
          direction_(dispersal_direction)
    {}

    /*! \copybrief RadialDispersalKernel::operator()()
     *
     * The randomness is based on the *generator*.
     * The new position does not depend on the position of the current
     * disperser thus *row* and *col* are unused.
     */
    template<typename Generator>
    std::tuple<int, int> operator() (Generator& generator, int row, int col)
    {
        switch (direction_) {
        case Direction::E:
            col += 1;
            break;
        case Direction::N:
            row -= 1;
            break;
        case Direction::NE:
            row -= 1;
            col += 1;
            break;
        case Direction::NW:
            row -= 1;
            col -= 1;
            break;
        case Direction::S:
            row += 1;
            break;
        case Direction::SE:
            row += 1;
            col += 1;
            break;
        case Direction::SW:
            row += 1;
            col -= 1;
            break;
        case Direction::W:
            col -= 1;
            break;
        default:
            throw std::invalid_argument("NeighborDispersalKernel: Unsupported Direction");
        }
        return std::make_tuple(row, col);
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        return type == DispersalKernelType::Uniform;
    }
};

} // namespace pops

#endif // POPS_NEIGHBOR_KERNEL_HPP
