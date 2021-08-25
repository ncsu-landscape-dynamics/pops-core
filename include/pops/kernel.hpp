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

#ifndef POPS_KERNEL_HPP
#define POPS_KERNEL_HPP

#include "switch_kernel.hpp"
#include "natural_anthropogenic_kernel.hpp"
#include "config.hpp"

namespace pops {

/**
 * @brief Create natural kernel
 *
 * Kernel parameters are taken from the configuration.
 *
 * @param dispersers The disperser raster (reference, for deterministic kernel)
 * @param network Network (initialized or not)
 * @return Created kernel
 */
template<typename IntegerRaster, typename RasterIndex>
SwitchDispersalKernel<IntegerRaster, RasterIndex> create_static_natural_kernel(
    const Config& config,
    const IntegerRaster& dispersers,
    const Network<RasterIndex>& network)
{
    auto natural_kernel = kernel_type_from_string(config.natural_kernel_type);
    UniformDispersalKernel uniform_kernel(config.rows, config.cols);
    RadialDispersalKernel<IntegerRaster> radial_kernel(
        config.ew_res,
        config.ns_res,
        natural_kernel,
        config.natural_scale,
        direction_from_string(config.natural_direction),
        config.natural_kappa,
        config.shape);
    DeterministicNeighborDispersalKernel natural_neighbor_kernel(
        direction_from_string(config.natural_direction));
    DeterministicDispersalKernel<IntegerRaster> deterministic_kernel(
        natural_kernel,
        dispersers,
        config.dispersal_percentage,
        config.ew_res,
        config.ns_res,
        config.natural_scale,
        config.shape);
    NetworkDispersalKernel<RasterIndex> network_kernel(
        network, config.network_min_time, config.network_max_time);
    SwitchDispersalKernel<IntegerRaster, RasterIndex> selectable_kernel(
        natural_kernel,
        radial_kernel,
        deterministic_kernel,
        uniform_kernel,
        network_kernel,
        natural_neighbor_kernel,
        config.deterministic);
    return selectable_kernel;
}
/**
 * @brief Create anthropogenic kernel
 *
 * Same structure as the natural kernel, but the parameters are for anthropogenic
 * kernel when available.
 *
 * @param dispersers The disperser raster (reference, for deterministic kernel)
 * @param network Network (initialized or not)
 * @return Created kernel
 */
template<typename IntegerRaster, typename RasterIndex>
SwitchDispersalKernel<IntegerRaster, RasterIndex> create_static_anthro_kernel(
    const Config& config,
    const IntegerRaster& dispersers,
    const Network<RasterIndex>& network)
{
    auto anthro_kernel = kernel_type_from_string(config.anthro_kernel_type);
    UniformDispersalKernel uniform_kernel(config.rows, config.cols);
    RadialDispersalKernel<IntegerRaster> radial_kernel(
        config.ew_res,
        config.ns_res,
        anthro_kernel,
        config.anthro_scale,
        direction_from_string(config.anthro_direction),
        config.anthro_kappa,
        config.shape);
    DeterministicNeighborDispersalKernel anthro_neighbor_kernel(
        direction_from_string(config.anthro_direction));
    DeterministicDispersalKernel<IntegerRaster> deterministic_kernel(
        anthro_kernel,
        dispersers,
        config.dispersal_percentage,
        config.ew_res,
        config.ns_res,
        config.anthro_scale,
        config.shape);
    NetworkDispersalKernel<RasterIndex> network_kernel(
        network, config.network_min_time, config.network_max_time);
    SwitchDispersalKernel<IntegerRaster, RasterIndex> selectable_kernel(
        anthro_kernel,
        radial_kernel,
        deterministic_kernel,
        uniform_kernel,
        network_kernel,
        anthro_neighbor_kernel,
        config.deterministic);
    return selectable_kernel;
}

// TODO: Use smart pointers or ensure delete is called.

/**
 * @brief Create natural kernel
 *
 * Kernel parameters are taken from the configuration.
 *
 * @param dispersers The disperser raster (reference, for deterministic kernel)
 * @param network Network (initialized or not)
 * @return Created kernel
 */
template<typename Generator, typename IntegerRaster, typename RasterIndex>
KernelInterface<Generator>* create_dynamic_natural_kernel(
    const Config& config,
    const IntegerRaster& dispersers,
    const Network<RasterIndex>& network)
{
    auto natural_kernel = kernel_type_from_string(config.natural_kernel_type);
    if (natural_kernel == DispersalKernelType::Uniform) {
        UniformDispersalKernel kernel(config.rows, config.cols);
        return new AlwaysEligibleDynamicKernel<UniformDispersalKernel, Generator>(
            kernel);
    }
    else if (natural_kernel == DispersalKernelType::DeterministicNeighbor) {
        DeterministicNeighborDispersalKernel kernel(
            direction_from_string(config.natural_direction));
        return new AlwaysEligibleDynamicKernel<
            DeterministicNeighborDispersalKernel,
            Generator>(kernel);
    }
    else if (natural_kernel == DispersalKernelType::Network) {
        NetworkDispersalKernel<RasterIndex> network_kernel(
            network, config.network_min_time, config.network_max_time);
        return new DynamicKernel<NetworkDispersalKernel<RasterIndex>, Generator>(
            network_kernel);
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
        return new AlwaysEligibleDynamicKernel<
            DeterministicDispersalKernel<IntegerRaster>,
            Generator>(deterministic_kernel);
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
        return new AlwaysEligibleDynamicKernel<
            RadialDispersalKernel<IntegerRaster>,
            Generator>(radial_kernel);
    }
}

/**
 * @brief Create anthropogenic kernel
 *
 * Same structure as the natural kernel, but the parameters are for anthropogenic
 * kernel when available.
 *
 * @param dispersers The disperser raster (reference, for deterministic kernel)
 * @param network Network (initialized or not)
 * @return Created kernel
 */
template<typename Generator, typename IntegerRaster, typename RasterIndex>
KernelInterface<Generator>* create_dynamic_anthro_kernel(
    const Config& config,
    const IntegerRaster& dispersers,
    const Network<RasterIndex>& network)
{
    auto anthro_kernel = kernel_type_from_string(config.anthro_kernel_type);
    if (anthro_kernel == DispersalKernelType::Uniform) {
        UniformDispersalKernel kernel(config.rows, config.cols);
        return new AlwaysEligibleDynamicKernel<UniformDispersalKernel, Generator>(
            kernel);
    }
    else if (anthro_kernel == DispersalKernelType::DeterministicNeighbor) {
        DeterministicNeighborDispersalKernel kernel(
            direction_from_string(config.anthro_direction));
        return new AlwaysEligibleDynamicKernel<
            DeterministicNeighborDispersalKernel,
            Generator>(kernel);
    }
    else if (anthro_kernel == DispersalKernelType::Network) {
        NetworkDispersalKernel<RasterIndex> network_kernel(
            network, config.network_min_time, config.network_max_time);
        return new DynamicKernel<NetworkDispersalKernel<RasterIndex>, Generator>(
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
        return new AlwaysEligibleDynamicKernel<
            DeterministicDispersalKernel<IntegerRaster>,
            Generator>(deterministic_kernel);
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
        return new AlwaysEligibleDynamicKernel<
            RadialDispersalKernel<IntegerRaster>,
            Generator>(radial_kernel);
    }
}

/*! Dispersal kernel supporting all available kernels for natural and anthropogenic
 * distance spread.
 *
 * This is a typedef defining the main disperal kernel class in the PoPS
 * library. When you are using the library, this is default choice.
 *
 * The choices which are allowed by this class are:
 * ```
 *                     kernel type
 *                    /           \
 *                   /             \
 *              natural        anthropogenic
 *             /   |   \          /   |   \
 *        radial uniform ...   radial uniform ...
 *        /  | \_____                 |
 *       /   |       \                |
 * cauchy exponential ...            ...
 * ```
 * The ellipses represent whatever the SwitchDispersalKernel and
 * RadialDispersalKernel classes support.
 *
 * See NaturalAnthropogenicDispersalKernel and SwitchDispersalKernel for further
 * documentation.
 */
template<typename IntegerRaster, typename RasterIndex>
using DispersalKernel = NaturalAnthropogenicDispersalKernel<
    SwitchDispersalKernel<IntegerRaster, RasterIndex>,
    SwitchDispersalKernel<IntegerRaster, RasterIndex>>;

template<typename Generator>
using DynamicDispersalKernel = DynamicNaturalAnthropogenicDispersalKernel<
    KernelInterface<Generator>,
    KernelInterface<Generator>>;

template<typename IntegerRaster, typename RasterIndex>
DispersalKernel<IntegerRaster, RasterIndex> create_static_kernel(
    const Config& config,
    const IntegerRaster& dispersers,
    const Network<RasterIndex>& network)
{
    return DispersalKernel<IntegerRaster, RasterIndex>(
        create_static_natural_kernel(config, dispersers, network),
        create_static_anthro_kernel(config, dispersers, network),
        config.use_anthropogenic_kernel,
        config.percent_natural_dispersal);
}

template<typename Generator, typename IntegerRaster, typename RasterIndex>
DynamicDispersalKernel<std::default_random_engine> create_dynamic_kernel(
    const Config& config,
    const IntegerRaster& dispersers,
    const Network<RasterIndex>& network)
{
    return DynamicDispersalKernel<Generator>(
        create_dynamic_natural_kernel<Generator>(config, dispersers, network),
        create_dynamic_anthro_kernel<Generator>(config, dispersers, network),
        config.use_anthropogenic_kernel,
        config.percent_natural_dispersal);
}

}  // namespace pops

#endif  // POPS_KERNEL_HPP
