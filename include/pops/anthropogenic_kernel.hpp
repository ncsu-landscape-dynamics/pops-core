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

#ifndef POPS_ANTHROPOGENIC_KERNEL_HPP
#define POPS_ANTHROPOGENIC_KERNEL_HPP

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
 * @brief Create anthropogenic kernel from configuration
 *
 * Same structure as the natural kernel, but the parameters are for anthropogenic
 * kernel when available.
 *
 * @param dispersers The disperser raster (reference, for deterministic kernel)
 * @param network Network (initialized or not)
 * @return Created kernel
 */
template<typename Generator, typename IntegerRaster, typename RasterIndex>
std::unique_ptr<KernelInterface<Generator>> create_anthro_kernel(
    const Config& config,
    const IntegerRaster& dispersers,
    const Network<RasterIndex>& network)
{
    auto anthro_kernel = kernel_type_from_string(config.anthro_kernel_type);
    if (anthro_kernel == DispersalKernelType::Uniform) {
        UniformDispersalKernel kernel(config.rows, config.cols);
        return std::make_unique<
            AlwaysEligibleDynamicKernel<UniformDispersalKernel, Generator>>(kernel);
    }
    else if (anthro_kernel == DispersalKernelType::DeterministicNeighbor) {
        DeterministicNeighborDispersalKernel kernel(
            direction_from_string(config.anthro_direction));
        return std::make_unique<AlwaysEligibleDynamicKernel<
            DeterministicNeighborDispersalKernel,
            Generator>>(kernel);
    }
    else if (anthro_kernel == DispersalKernelType::Network) {
        NetworkDispersalKernel<RasterIndex> network_kernel(
            network, config.network_min_time, config.network_max_time);
        return std::make_unique<
            DynamicKernel<NetworkDispersalKernel<RasterIndex>, Generator>>(
            network_kernel);
    }
    else if (config.deterministic) {
        DeterministicDispersalKernel<IntegerRaster> deterministic_kernel(
            anthro_kernel,
            dispersers,
            config.dispersal_percentage,
            config.ew_res,
            config.ns_res,
            config.anthro_scale,
            config.shape);
        return std::make_unique<AlwaysEligibleDynamicKernel<
            DeterministicDispersalKernel<IntegerRaster>,
            Generator>>(deterministic_kernel);
    }
    else {
        RadialDispersalKernel<IntegerRaster> radial_kernel(
            config.ew_res,
            config.ns_res,
            anthro_kernel,
            config.anthro_scale,
            direction_from_string(config.anthro_direction),
            config.anthro_kappa,
            config.shape);
        return std::make_unique<AlwaysEligibleDynamicKernel<
            RadialDispersalKernel<IntegerRaster>,
            Generator>>(radial_kernel);
    }
}

}  // namespace pops

#endif  // POPS_ANTHROPOGENIC_KERNEL_HPP
