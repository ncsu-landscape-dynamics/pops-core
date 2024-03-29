/*
 * PoPS model - dispersal kernel combining natural and anthropogenic dispersals
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

#ifndef POPS_NATURAL_ANTHROPOGENIC_KERNEL_HPP
#define POPS_NATURAL_ANTHROPOGENIC_KERNEL_HPP

#include "kernel_types.hpp"

#include <memory>
#include <tuple>
#include <random>
#include <type_traits>

namespace pops {

/*! Dispersal kernel template for natural and anthropogenic distance dispersal
 *
 * This class template can be used for switching between two kernels
 * specified as template parameters where one kernel is assumed to
 * be used natural distance distance dispersal and the other for anthropogenic
 * distance one. The two dispersal kernels can be the same or different.
 *
 * There are no assumptions for natural or anthropogenic kernels. We are using
 * these names for clarity since this is how we expect the template
 * to be used.
 *
 * Bernoulli distribution is used to decide between the natural and
 * anthropogenic distance kernel. The anthropogenic distance dispersal can be also
 * competely disabled.
 *
 * Instances cannot be copied due to contained std::unique_ptr, but they can be
 * moved, so returning new objects from functions works. This is enforced by the
 * compiler.
 */

template<typename NaturalKernelType, typename AnthropogenicKernelType>
class NaturalAnthropogenicDispersalKernel
{
protected:
    bool use_anthropogenic_kernel_;
    std::unique_ptr<NaturalKernelType> natural_kernel_;
    std::unique_ptr<AnthropogenicKernelType> anthropogenic_kernel_;
    std::bernoulli_distribution bernoulli_distribution;

public:
    NaturalAnthropogenicDispersalKernel(
        std::unique_ptr<NaturalKernelType> natural_kernel,
        std::unique_ptr<AnthropogenicKernelType> anthropogenic_kernel,
        bool use_anthropogenic_kernel,
        double percent_natural_dispersal)
        : use_anthropogenic_kernel_(use_anthropogenic_kernel),
          // Here we initialize all distributions,
          // although we won't use all of them.
          natural_kernel_(std::move(natural_kernel)),
          anthropogenic_kernel_(std::move(anthropogenic_kernel)),
          // use bernoulli distribution to act as the sampling with prob(gamma,1-gamma)
          bernoulli_distribution(percent_natural_dispersal)
    {}

    /*! \copydoc RadialDispersalKernel::operator()()
     */
    template<typename Generator>
    std::tuple<int, int> operator()(Generator& generator, int row, int col)
    {
        // switch in between the supported kernels
        if (!use_anthropogenic_kernel_
            || !anthropogenic_kernel_->is_cell_eligible(row, col)
            || bernoulli_distribution(generator.anthropogenic_dispersal())) {
            return natural_kernel_->operator()(generator.natural_dispersal(), row, col);
        }
        return anthropogenic_kernel_->operator()(
            generator.anthropogenic_dispersal(), row, col);
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     *
     * Returns true if at least one of the kernels (natural or anthropogenic)
     * supports the given kernel type.
     *
     * \note Note that if natural and anthropogenic kernels are different, this is
     * not generally usable because one kernel can support that and the
     * other not. However, there is not much room for accidental misuse
     * of this because this class does not use the type directly
     * (it is handled by the underlying kernels).
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        if (std::is_same<NaturalKernelType, AnthropogenicKernelType>::value) {
            return NaturalKernelType::supports_kernel(type);
        }
        else {
            return NaturalKernelType::supports_kernel(type)
                   || AnthropogenicKernelType::supports_kernel(type);
        }
    }
};

}  // namespace pops

#endif  // POPS_NATURAL_ANTHROPOGENIC_KERNEL_HPP
