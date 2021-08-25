/*
 * PoPS model - disperal kernels
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

#ifndef POPS_SWITCH_KERNEL_HPP
#define POPS_SWITCH_KERNEL_HPP

#include "radial_kernel.hpp"
#include "deterministic_kernel.hpp"
#include "uniform_kernel.hpp"
#include "neighbor_kernel.hpp"
#include "network_kernel.hpp"
#include "kernel_types.hpp"
#include "utils.hpp"

#include <any>

namespace pops {

template<typename Generator>
class KernelInterface
{
public:
    /*! \copydoc RadialDispersalKernel::operator()()
     */
    virtual std::tuple<int, int> operator()(Generator& generator, int row, int col) = 0;

    virtual bool is_cell_eligible(int row, int col) = 0;

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    virtual bool supports_kernel(const DispersalKernelType type) = 0;

    virtual ~KernelInterface() = default;
};

template<typename ActualKernel, typename Generator>
class BaseDynamicKernel : public KernelInterface<Generator>
{
public:
    BaseDynamicKernel(const ActualKernel& kernel) : kernel_(kernel) {}

    /*! \copydoc RadialDispersalKernel::operator()()
     */
    std::tuple<int, int> operator()(Generator& generator, int row, int col) override
    {
        return kernel_.operator()(generator, row, col);
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    bool supports_kernel(const DispersalKernelType type) override
    {
        return ActualKernel::supports_kernel(type);
    }

protected:
    ActualKernel kernel_;
};

template<typename ActualKernel, typename Generator>
class DynamicKernel : public BaseDynamicKernel<ActualKernel, Generator>
{
public:
    DynamicKernel(const ActualKernel& kernel)
        : BaseDynamicKernel<ActualKernel, Generator>(kernel)
    {}

    bool is_cell_eligible(int row, int col) override
    {
        return BaseDynamicKernel<ActualKernel, Generator>::kernel_.is_cell_eligible(
            row, col);
    }
};

template<typename ActualKernel, typename Generator>
class AlwaysEligibleDynamicKernel : public BaseDynamicKernel<ActualKernel, Generator>
{
public:
    AlwaysEligibleDynamicKernel(const ActualKernel& kernel)
        : BaseDynamicKernel<ActualKernel, Generator>(kernel)
    {}

    bool is_cell_eligible(int row, int col) override
    {
        UNUSED(row);
        UNUSED(col);
        return true;
    }
};

/*! Dispersal kernel providing all the radial kernels.
 *
 * We understand a radial kernel to be a kernel which has parameters
 * which translate into a distance and direction.
 *
 * To add new kernel, add new member, constructor parameter,
 * its call in the function call operator, and extend the
 * supports_kernel() function.
 */
// template<typename IntegerRaster, typename RasterIndex>
// TODO: unique pointer here
/*
KernelInterface* create_dynamic_kernel(
    const DispersalKernelType& dispersal_kernel_type, const bool deterministic = false)
{
    if (dispersal_kernel_type_ == DispersalKernelType::Uniform) {
        return DynamicKernel(UniformDispersalKernel());
    }
    else if (dispersal_kernel_type_ == DispersalKernelType::DeterministicNeighbor) {
        return DeterministicNeighborDispersalKernel(Direction::None);
    }
    else if (dispersal_kernel_type_ == DispersalKernelType::Network) {
        return NetworkDispersalKernel<RasterIndex>();
    }
    else if (deterministic_) {
        return DeterministicDispersalKernel<IntegerRaster>();
    }
    else {
        return RadialDispersalKernel<IntegerRaster>(generator, row, col);
    }
}
*/

/*! Dispersal kernel providing all the radial kernels.
 *
 * We understand a radial kernel to be a kernel which has parameters
 * which translate into a distance and direction.
 *
 * To add new kernel, add new member, constructor parameter,
 * its call in the function call operator, and extend the
 * supports_kernel() function.
 */
template<typename IntegerRaster, typename RasterIndex>
class SwitchDispersalKernel
{
protected:
    DispersalKernelType dispersal_kernel_type_;
    RadialDispersalKernel<IntegerRaster> radial_kernel_;
    DeterministicDispersalKernel<IntegerRaster> deterministic_kernel_;
    UniformDispersalKernel uniform_kernel_;
    DeterministicNeighborDispersalKernel deterministic_neighbor_kernel_;
    NetworkDispersalKernel<RasterIndex> network_kernel_;
    bool deterministic_;

public:
    SwitchDispersalKernel(
        const DispersalKernelType& dispersal_kernel_type,
        const RadialDispersalKernel<IntegerRaster>& radial_kernel,
        const DeterministicDispersalKernel<IntegerRaster>& deterministic_kernel,
        const UniformDispersalKernel& uniform_kernel,
        const NetworkDispersalKernel<RasterIndex>& network_kernel,
        const DeterministicNeighborDispersalKernel& deterministic_neighbor_kernel =
            DeterministicNeighborDispersalKernel(Direction::None),
        const bool deterministic = false)
        : dispersal_kernel_type_(dispersal_kernel_type),
          // Here we initialize all kernels,
          // although we won't use all of them.
          radial_kernel_(radial_kernel),
          deterministic_kernel_(deterministic_kernel),
          uniform_kernel_(uniform_kernel),
          deterministic_neighbor_kernel_(deterministic_neighbor_kernel),
          network_kernel_(network_kernel),
          deterministic_(deterministic)
    {}

    /*! \copydoc RadialDispersalKernel::operator()()
     */
    template<typename Generator>
    std::tuple<int, int> operator()(Generator& generator, int row, int col)
    {
        // switch in between the supported kernels
        if (dispersal_kernel_type_ == DispersalKernelType::Uniform) {
            return uniform_kernel_(generator, row, col);
        }
        else if (dispersal_kernel_type_ == DispersalKernelType::DeterministicNeighbor) {
            return deterministic_neighbor_kernel_(generator, row, col);
        }
        else if (dispersal_kernel_type_ == DispersalKernelType::Network) {
            return network_kernel_(generator, row, col);
        }
        else if (deterministic_) {
            return deterministic_kernel_(generator, row, col);
        }
        else {
            return radial_kernel_(generator, row, col);
        }
    }

    bool is_cell_eligible(int row, int col)
    {
        // switch in between the supported kernels
        if (dispersal_kernel_type_ == DispersalKernelType::Uniform) {
            // TODO: Individual kernels should support this.
            return true;
        }
        else if (dispersal_kernel_type_ == DispersalKernelType::DeterministicNeighbor) {
            return true;
        }
        else if (dispersal_kernel_type_ == DispersalKernelType::Network) {
            return network_kernel_.is_cell_eligible(row, col);
        }
        else if (deterministic_) {
            return true;
        }
        else {
            return true;
        }
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        if (type == DispersalKernelType::Uniform) {
            return true;
        }
        else if (type == DispersalKernelType::DeterministicNeighbor) {
            return true;
        }
        // Radial and Deterministic support the same kernel types
        else {
            return RadialDispersalKernel<IntegerRaster>::supports_kernel(type);
        }
    }
};

}  // namespace pops

#endif  // POPS_SWITCH_KERNEL_HPP
