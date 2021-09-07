/*
 * PoPS model - disperal kernels
 *
 * Copyright (C) 2021 by the authors.
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

#ifndef KERNEL_BASE_HPP
#define KERNEL_BASE_HPP

#include "kernel_types.hpp"

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

}  // namespace pops

#endif  // KERNEL_BASE_HPP
