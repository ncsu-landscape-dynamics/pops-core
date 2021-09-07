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

#ifndef POPS_NATURAL_KERNEL_HPP
#define POPS_NATURAL_KERNEL_HPP

#include "radial_kernel.hpp"
#include "deterministic_kernel.hpp"
#include "uniform_kernel.hpp"
#include "neighbor_kernel.hpp"
#include "kernel_types.hpp"
#include "kernel_base.hpp"
#include "config.hpp"

#include <memory>

namespace pops {

/**
 * @brief Create natural kernel from configuration
 *
 * Kernel parameters are taken from the configuration.
 *
 * @param dispersers The disperser raster (reference, for deterministic kernel)
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
