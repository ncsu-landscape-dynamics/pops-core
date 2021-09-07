/*
 * PoPS model - main disperal kernel
 *
 * Copyright (C) 2019-2021 by the authors.
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

/*! \file kernel.hpp
 *
 * \brief Main entry point to dispersal kernel functionality
 *
 * This file contains convenient definitions or wrappers to be used
 * when using this library.
 *
 * ## Adding a new kernel
 *
 * To add new kernel, decide if it needs to be a separate class,
 * or if it is just a parameterization of the existing kernel.
 * Many kernels can be handled by the RadialDispersalKernel class.
 *
 * Generally, such class needs to provide contructor taking all required
 * parameters and a function call operator which takes a random number
 * generator object and returns a pair of coordinates as row and column.
 * The signature of the function template should be:
 *
 * ```
 * template<typename Generator>
 * std::tuple<int, int> operator() (Generator& generator)
 * ```
 *
 * Besides implementation in a class, enum for the different types
 * of kernels needs to be extented as well as function which transforms
 * strings into enum values, i.e., you need to add to the
 * `DispersalKernelType` enum and extend the kernel_type_from_string()
 * function in kernel_types.hpp.
 */

#ifndef POPS_NATURAL_KERNEL_HPP
#define POPS_NATURAL_KERNEL_HPP

#include "radial_kernel.hpp"
#include "deterministic_kernel.hpp"
#include "uniform_kernel.hpp"
#include "neighbor_kernel.hpp"
#include "network_kernel.hpp"
#include "kernel_types.hpp"
#include "kernel_base.hpp"
#include "natural_anthropogenic_kernel.hpp"
#include "config.hpp"

#include <memory>

namespace pops {

/**
 * @brief Create natural kernel from configuration
 *
 * Kernel parameters are taken from the configuration.
 *
 * @param dispersers The disperser raster (reference, for deterministic kernel)
 * @param network Network (initialized or not)
 * @return Created kernel
 */
template<typename Generator, typename IntegerRaster, typename RasterIndex>
std::unique_ptr<KernelInterface<Generator>>
create_natural_kernel(const Config& config, const IntegerRaster& dispersers)
{
    auto natural_kernel = kernel_type_from_string(config.natural_kernel_type);
    if (natural_kernel == DispersalKernelType::Uniform) {
        UniformDispersalKernel kernel(config.rows, config.cols);
        return std::make_unique<
            AlwaysEligibleDynamicKernel<UniformDispersalKernel, Generator>>(kernel);
    }
    else if (natural_kernel == DispersalKernelType::DeterministicNeighbor) {
        DeterministicNeighborDispersalKernel kernel(
            direction_from_string(config.natural_direction));
        return std::make_unique<AlwaysEligibleDynamicKernel<
            DeterministicNeighborDispersalKernel,
            Generator>>(kernel);
    }
    else if (config.deterministic) {

        DeterministicDispersalKernel<IntegerRaster> deterministic_kernel(
            natural_kernel,
            dispersers,
            config.dispersal_percentage,
            config.ew_res,
            config.ns_res,
            config.natural_scale,
            config.shape);
        return std::make_unique<AlwaysEligibleDynamicKernel<
            DeterministicDispersalKernel<IntegerRaster>,
            Generator>>(deterministic_kernel);
    }
    else {
        RadialDispersalKernel<IntegerRaster> radial_kernel(
            config.ew_res,
            config.ns_res,
            natural_kernel,
            config.natural_scale,
            direction_from_string(config.natural_direction),
            config.natural_kappa,
            config.shape);
        return std::make_unique<AlwaysEligibleDynamicKernel<
            RadialDispersalKernel<IntegerRaster>,
            Generator>>(radial_kernel);
    }
}

}  // namespace pops

#endif  // POPS_NATURAL_KERNEL_HPP
