#ifdef POPS_TEST

/*
 * Test of PoPS kernels with simulation class focused on benchmarking.
 *
 * Copyright (C) 2021 by the authors.
 *
 * Authors: Vaclav Petras <wenzeslaus gmail com>
 *
 * This file is part of PoPS.

 * PoPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * PoPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with PoPS. If not, see <https://www.gnu.org/licenses/>.
 */

#include <pops/config.hpp>
#include <pops/raster.hpp>
#include <pops/utils.hpp>
#include <pops/kernel.hpp>
#include <pops/switch_kernel.hpp>
#include <pops/simulation.hpp>
#include <pops/model.hpp>

#include <map>
#include <iostream>
#include <fstream>
#include <vector>

#include <regex>

using namespace pops;

template<typename NaturalKernelType, typename AnthropogenicKernelType>
class StaticNaturalAnthropogenicDispersalKernel
{
protected:
    bool use_anthropogenic_kernel_;
    NaturalKernelType natural_kernel_;
    AnthropogenicKernelType anthropogenic_kernel_;
    std::bernoulli_distribution bernoulli_distribution;

public:
    StaticNaturalAnthropogenicDispersalKernel(
        const NaturalKernelType& natural_kernel,
        const AnthropogenicKernelType& anthropogenic_kernel,
        bool use_anthropogenic_kernel,
        double percent_natural_dispersal)
        : use_anthropogenic_kernel_(use_anthropogenic_kernel),
          // Here we initialize all distributions,
          // although we won't use all of them.
          natural_kernel_(natural_kernel),
          anthropogenic_kernel_(anthropogenic_kernel),
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
            || !anthropogenic_kernel_.is_cell_eligible(row, col)
            || bernoulli_distribution(generator.anthropogenic_dispersal())) {
            return natural_kernel_(generator.natural_dispersal(), row, col);
        }
        return anthropogenic_kernel_(generator.anthropogenic_dispersal(), row, col);
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

/**
 * @brief Create natural kernel from configuration
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
        network, config.network_min_distance, config.network_max_distance);
    SwitchDispersalKernel<IntegerRaster, RasterIndex> selectable_kernel(
        natural_kernel,
        radial_kernel,
        deterministic_kernel,
        uniform_kernel,
        network_kernel,
        natural_neighbor_kernel,
        config.dispersal_stochasticity);
    return selectable_kernel;
}
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
        network, config.network_min_distance, config.network_max_distance);
    SwitchDispersalKernel<IntegerRaster, RasterIndex> selectable_kernel(
        anthro_kernel,
        radial_kernel,
        deterministic_kernel,
        uniform_kernel,
        network_kernel,
        anthro_neighbor_kernel,
        config.dispersal_stochasticity);
    return selectable_kernel;
}

template<typename IntegerRaster, typename RasterIndex>
using StaticDispersalKernel = StaticNaturalAnthropogenicDispersalKernel<
    SwitchDispersalKernel<IntegerRaster, RasterIndex>,
    SwitchDispersalKernel<IntegerRaster, RasterIndex>>;

template<typename IntegerRaster, typename RasterIndex>
StaticDispersalKernel<IntegerRaster, RasterIndex> create_static_kernel(
    const Config& config,
    const IntegerRaster& dispersers,
    const Network<RasterIndex>& network)
{
    return StaticDispersalKernel<IntegerRaster, RasterIndex>(
        create_static_natural_kernel(config, dispersers, network),
        create_static_anthro_kernel(config, dispersers, network),
        config.use_anthropogenic_kernel,
        config.percent_natural_dispersal);
}

template<typename KernelFactory>
int test_simulation_with_kernels_generic(
    const KernelFactory& kernel_factory, int steps, const std::string& kernel_type)
{
    int rows = 500;
    int cols = 1000;

    Config config;
    config.rows = rows;
    config.cols = cols;
    config.use_treatments = false;
    config.ew_res = 30;
    config.ns_res = 30;

    config.weather = false;
    config.reproductive_rate = 2;
    config.generate_stochasticity = false;
    config.establishment_stochasticity = false;
    std::vector<std::vector<int>> movements;
    // We want everything to establish.
    config.establishment_probability = 1;
    config.natural_kernel_type = kernel_type;
    config.natural_direction = "none";
    config.natural_scale = 1;
    config.anthro_scale = 1;
    config.natural_kappa = 0;
    config.anthro_kappa = 0;
    config.dispersal_percentage = 0.99;
    config.use_anthropogenic_kernel = false;
    config.random_seed = 42;
    config.model_type = "SI";
    config.latency_period_steps = 3;
    config.use_lethal_temperature = false;
    config.use_survival_rate = false;
    config.use_quarantine = false;
    config.use_spreadrates = false;

    config.set_date_start(2020, 1, 1);
    config.set_date_end(2021, 12, 31);
    config.set_step_unit(StepUnit::Month);
    config.set_step_num_units(1);
    config.use_mortality = false;
    config.mortality_frequency = "year";
    config.mortality_frequency_n = 1;

    config.create_schedules();

    config.dispersal_stochasticity = true;

    using IntRaster = Raster<int, int>;
    using DoubleRaster = Raster<double, int>;
    IntRaster infected{rows, cols, 10};
    IntRaster mortality_tracker{rows, cols, 0};
    IntRaster total_exposed{rows, cols, 0};
    // Susceptible and total are set in a way that there won't be any
    // dilution effect and the disperser will always establish given the
    // selected random seed. Establishment probability is high and with
    // the given seed we don't get any random numbers in establishment
    // test higher than that. (The weather is disabled.)
    IntRaster susceptible{rows, cols, 100};
    // add a lot of hosts, so that exposing or infecting them won't
    // chanage the susceptible/total ratio much
    // we want to minimize the dilution effect
    IntRaster total_hosts = susceptible;
    std::vector<std::vector<int>> suitable_cells =
        find_suitable_cells<int>(susceptible);

    IntRaster dispersers(infected.rows(), infected.cols());
    IntRaster established_dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    std::vector<IntRaster> exposed(
        config.latency_period_steps + 1,
        IntRaster(infected.rows(), infected.cols(), 0));

    DefaultSingleGeneratorProvider generator(config.random_seed);
    Simulation<IntRaster, DoubleRaster, int> simulation(
        config.rows,
        config.cols,
        model_type_from_string(config.model_type),
        config.latency_period_steps);

    for (int i = 0; i < steps; ++i) {
        dispersers = config.reproductive_rate * infected;
        auto kernel = kernel_factory(config, dispersers, Network<int>::null_network());

        simulation.disperse_and_infect(
            i,
            dispersers,
            established_dispersers,
            susceptible,
            exposed,
            infected,
            mortality_tracker,
            total_hosts,
            total_exposed,
            outside_dispersers,
            config.weather,
            kernel,
            suitable_cells,
            0.5,
            generator);
    }
    return 0;
}

template<typename KernelFactory>
int test_model_with_kernels_generic(
    const KernelFactory& kernel_factory, int steps, const std::string& kernel_type)
{
    int rows = 500;
    int cols = 1000;

    Config config;
    config.rows = rows;
    config.cols = cols;
    config.use_treatments = false;
    config.ew_res = 30;
    config.ns_res = 30;

    config.weather = false;
    config.reproductive_rate = 2;
    config.generate_stochasticity = false;
    config.establishment_stochasticity = false;
    std::vector<std::vector<int>> movements;
    // We want everything to establish.
    config.establishment_probability = 1;
    config.natural_kernel_type = kernel_type;
    config.natural_direction = "none";
    config.natural_scale = 1;
    config.anthro_scale = 1;
    config.natural_kappa = 0;
    config.anthro_kappa = 0;
    config.dispersal_percentage = 0.99;
    config.use_anthropogenic_kernel = false;
    config.random_seed = 42;
    config.model_type = "SI";
    config.latency_period_steps = 3;
    config.use_lethal_temperature = false;
    config.use_survival_rate = false;
    config.use_quarantine = false;
    config.use_spreadrates = false;

    Date date{2020, 1, 1};
    config.set_date_start(date);
    date.add_days(steps - 1);
    config.set_date_end(date);
    config.set_step_unit(StepUnit::Day);

    config.set_step_num_units(1);
    config.use_mortality = false;
    config.mortality_frequency = "year";
    config.mortality_frequency_n = 1;

    config.create_schedules();

    config.dispersal_stochasticity = true;

    Raster<int> infected{rows, cols, 10};
    // Susceptible and total are set in a way that there won't be any
    // dilution effect and the disperser will always establish given the
    // selected random seed. Establishment probability is high and with
    // the given seed we don't get any random numbers in establishment
    // test higher than that. (The weather is disabled.)
    Raster<int> susceptible{rows, cols, 100};
    Raster<int> total_hosts = susceptible + infected;
    Raster<int> total_populations = total_hosts;
    Raster<int> zeros(infected.rows(), infected.cols(), 0);
    Raster<int> dispersers(infected.rows(), infected.cols());
    Raster<int> established_dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    std::vector<std::vector<int>> suitable_cells =
        find_suitable_cells<int>(susceptible);

    unsigned num_mortality_steps = 1;
    std::vector<Raster<int>> mortality_tracker(
        num_mortality_steps, Raster<int>(infected.rows(), infected.cols(), 0));

    Raster<int> died(infected.rows(), infected.cols(), 0);

    // Exposed
    Raster<int> total_exposed(infected.rows(), infected.cols(), 0);

    std::vector<Raster<int>> empty_integer;
    std::vector<Raster<double>> empty_float;
    Treatments<Raster<int>, Raster<double>> treatments(config.scheduler());
    QuarantineEscapeAction<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, 0);

    Model<
        Raster<int>,
        Raster<double>,
        Raster<double>::IndexType,
        std::default_random_engine,
        KernelFactory>
        model(config, kernel_factory);
    for (unsigned step = 0; step < config.scheduler().get_num_steps(); ++step) {
        model.run_step(
            step,
            infected,
            susceptible,
            total_populations,
            total_hosts,
            dispersers,
            established_dispersers,
            total_exposed,
            empty_integer,
            mortality_tracker,
            died,
            empty_float,
            empty_float,
            treatments,
            zeros,
            outside_dispersers,
            quarantine,
            zeros,
            movements,
            Network<int>::null_network(),
            suitable_cells);
    }
    return 0;
}

using SimpleDDKernel = StaticNaturalAnthropogenicDispersalKernel<
    DeterministicNeighborDispersalKernel,
    DeterministicNeighborDispersalKernel>;

template<typename IntegerRaster, typename RasterIndex>
SimpleDDKernel create_simple_static_kernel(
    const Config& config,
    const IntegerRaster& dispersers,
    const Network<RasterIndex>& network)
{
    UNUSED(config);
    UNUSED(dispersers);
    UNUSED(network);
    DeterministicNeighborDispersalKernel kernel(Direction::E);
    return SimpleDDKernel(kernel, kernel, false, 1);
}

template<typename IntegerRaster, typename RasterIndex>
RadialDispersalKernel<IntegerRaster> create_radial_kernel(
    const Config& config,
    const IntegerRaster& dispersers,
    const Network<RasterIndex>& network)
{
    UNUSED(dispersers);
    UNUSED(network);
    return RadialDispersalKernel<IntegerRaster>(
        config.ew_res,
        config.ns_res,
        kernel_type_from_string(config.natural_kernel_type),
        config.natural_scale,
        direction_from_string(config.natural_direction),
        config.natural_kappa,
        config.shape);
}

template<typename IntegerRaster>
using DoubleRadialKernel = StaticNaturalAnthropogenicDispersalKernel<
    RadialDispersalKernel<IntegerRaster>,
    RadialDispersalKernel<IntegerRaster>>;
template<typename IntegerRaster, typename RasterIndex>
DoubleRadialKernel<IntegerRaster> create_double_radial_kernel(
    const Config& config,
    const IntegerRaster& dispersers,
    const Network<RasterIndex>& network)
{
    return DoubleRadialKernel<IntegerRaster>(
        create_radial_kernel(config, dispersers, network),
        RadialDispersalKernel<IntegerRaster>(
            config.ew_res,
            config.ns_res,
            kernel_type_from_string(config.anthro_kernel_type),
            config.anthro_scale,
            direction_from_string(config.anthro_direction),
            config.anthro_kappa,
            config.shape),
        config.use_anthropogenic_kernel,
        config.percent_natural_dispersal);
}

int test_simulation_with_kernels(
    std::string command, int steps, const std::string& kernel_type)
{
    int ret = 0;
    if (command == "trivial") {
        ret += test_simulation_with_kernels_generic(
            create_simple_static_kernel<Raster<int>, int>, steps, kernel_type);
    }
    else if (command == "radial") {
        ret += test_simulation_with_kernels_generic(
            create_radial_kernel<Raster<int>, int>, steps, kernel_type);
    }
    else if (command == "double") {
        ret += test_simulation_with_kernels_generic(
            create_double_radial_kernel<Raster<int>, int>, steps, kernel_type);
    }
    else if (command == "static") {
        ret += test_simulation_with_kernels_generic(
            create_static_kernel<Raster<int>, int>, steps, kernel_type);
    }
    else if (command == "dynamic") {
        ret += test_simulation_with_kernels_generic(
            create_dynamic_kernel<std::default_random_engine, Raster<int>, int>,
            steps,
            kernel_type);
    }
    else {
        std::cerr << "Unknown sub-command: " << command << "\n";
        return 1;
    }
    return ret;
}

int test_model_with_kernels(
    std::string command, int steps, const std::string& kernel_type)
{
    int ret = 0;
    if (command == "trivial") {
        ret += test_model_with_kernels_generic(
            create_simple_static_kernel<Raster<int>, int>, steps, kernel_type);
    }
    else if (command == "radial") {
        ret += test_model_with_kernels_generic(
            create_radial_kernel<Raster<int>, int>, steps, kernel_type);
    }
    else if (command == "double") {
        ret += test_model_with_kernels_generic(
            create_double_radial_kernel<Raster<int>, int>, steps, kernel_type);
    }
    else if (command == "static") {
        ret += test_model_with_kernels_generic(
            create_static_kernel<Raster<int>, int>, steps, kernel_type);
    }
    else if (command == "dynamic") {
        ret += test_model_with_kernels_generic(
            create_dynamic_kernel<std::default_random_engine, Raster<int>, int>,
            steps,
            kernel_type);
    }
    else {
        std::cerr << "Unknown sub-command: " << command << "\n";
        return 1;
    }
    return ret;
}

int kernel_type_test(int argc, char** argv)
{
    if (argc != 5) {
        std::cerr
            << "Usage: " << argv[0]
            << "model|simulation trivial|radial|double|static|dynamic STEPS KERNEL_TYPE\n";
        return argc;
    }
    std::string command = argv[1];
    std::string subcommand = argv[2];
    int steps = std::stoi(argv[3]);
    std::string kernel_type = argv[4];

    int ret = 0;

    if (command == "model") {
        ret += test_model_with_kernels(subcommand, steps, kernel_type);
    }
    else if (command == "simulation") {
        ret += test_simulation_with_kernels(subcommand, steps, kernel_type);
    }
    else {
        std::cerr << "Unknown sub-command: " << command << "\n";
        return 1;
    }
    return ret;
}

int test_fixed_kernel_types()
{
    int ret = 0;
    int steps = 2;
    std::string kernel_type{"exponential"};
    for (std::string command : {"trivial", "dynamic"}) {
        ret += test_model_with_kernels(command, steps, kernel_type);
        ret += test_simulation_with_kernels(command, steps, kernel_type);
    }
    return ret;
}

int main(int argc, char** argv)
{
    int ret = 0;

    if (argc < 2)
        ret += test_fixed_kernel_types();
    else
        ret += kernel_type_test(argc, argv);

    return ret;
}

#endif  // POPS_TEST
