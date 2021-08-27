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
#include <pops/simulation.hpp>
#include <pops/model.hpp>

#include <map>
#include <iostream>
#include <fstream>
#include <vector>

#include <regex>

using namespace pops;

/** Placeholder for string to string no-op. */
std::string convert_to(const std::string& text, std::string tag)
{
    UNUSED(tag);
    return text;
}

/** Convert string to double */
double convert_to(const std::string& text, double tag)
{
    UNUSED(tag);
    return std::stod(text);
}

/** Convert string to int */
int convert_to(const std::string& text, int tag)
{
    UNUSED(tag);
    return std::stoi(text);
}

/**
 * \brief A generic key-value configuration which can convert values to desired types.
 */
class RawConfig
{
public:
    /**
     * Get value from config or the default value.
     *
     * The type is determined by the template parameter is specified.
     */
    template<typename T = std::string>
    T get(const std::string& key) const
    {
        auto it{values_.find(key)};
        if (it != values_.end())
            return convert_to(it->second, T());
        throw std::invalid_argument(std::string("No value for key: ") + key);
    }
    /**
     * Get value from config or the default value.
     *
     * The type is determined by the default value unless the template parameter is
     * specified. For floating point numbers, you can do:
     *
     * ```
     * double x = config.get("x", 0.);
     * ```
     *
     * However, if `0` is passed, the function parse the value as an integer and return
     * an integer.
     */
    template<typename T = std::string>
    T get(const std::string& key, T default_value) const
    {
        auto it{values_.find(key)};
        if (it != values_.end())
            return convert_to(it->second, T());
        return default_value;
    }
    // std::optional for default_value can replace the get overload in C++17.
    template<typename T>
    void set(const std::string& key, const T& value)
    {
        values_[key] = value;
    }

protected:
    std::map<std::string, std::string> values_;
};

/**
 * \brief Create RawConfig from a YAML-like text stream of simple `key: value` pairs.
 */
template<typename Stream>
RawConfig read_config(Stream& stream)
{
    RawConfig config;
    std::string line;
    while (std::getline(stream, line)) {
        std::regex delimeter(R"([\s]*:[\s]*)");
        std::smatch match;
        if (regex_search(line, match, delimeter)) {
            config.set(match.prefix(), match.suffix());
        }
        else {
            throw std::runtime_error(std::string("Incorrect format at line: ") + line);
        }
    }
    return config;
}

// void update_config_from_raw_config(Config& config, const RawConfig& raw) {}

template<typename Value>
void get_or_use(RawConfig raw, Value& val, const std::string& key)
{
    val = raw.get(key, val);
}

int test_kernels(int argc, char** argv)
{
    int rows = 500;
    int cols = 1000;
    using IntRaster = Raster<int, int>;
    using DoubleRaster = Raster<double, int>;
    IntRaster infected{rows, cols, 5};
    IntRaster mortality_tracker{rows, cols, 0};
    IntRaster total_exposed{rows, cols, 0};
    // Susceptible and total are set in a way that there won't be any
    // dilution effect and the disperser will always establish given the
    // selected random seed. Establishment probability is high and with
    // the given seed we don't get any random numbers in establishment
    // test higher than that. (The weather is disabled.)
    IntRaster susceptible{rows, cols, 1000};
    // add a lot of hosts, so that exposing or infecting them won't
    // chanage the susceptible/total ratio much
    // we want to minimize the dilution effect
    IntRaster total_hosts = susceptible;
    DoubleRaster temperature{rows, cols, 5};
    DoubleRaster weather_coefficient{rows, cols, 0};
    IntRaster zeros(infected.rows(), infected.cols(), 0);
    std::vector<std::vector<int>> suitable_cells =
        find_suitable_cells<int>(susceptible);

    IntRaster dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;
    bool weather = false;
    double reproductive_rate = 2;
    unsigned latency_period_steps = 3;

    std::vector<IntRaster> exposed(
        latency_period_steps + 1, IntRaster(infected.rows(), infected.cols(), 0));

    DeterministicNeighborDispersalKernel kernel(Direction::E);
    NaturalAnthropogenicDispersalKernel<
        DeterministicNeighborDispersalKernel,
        DeterministicNeighborDispersalKernel>
        kernel_s(kernel, kernel, false, 1);

    Simulation<IntRaster, DoubleRaster, int, std::default_random_engine> simulation(
        42,
        infected.rows(),
        infected.cols(),
        model_type_from_string("SEI"),
        latency_period_steps);

    /*
    DynamicDispersalKernel<std::default_random_engine> dispersal_kernel(
        create_natural_kernel(dispersers, network),
        create_anthro_kernel(dispersers, network),
        config_.use_anthropogenic_kernel,
        config_.percent_natural_dispersal);
    */

    KernelInterface<std::default_random_engine>* kernel2 =
        new AlwaysEligibleDynamicKernel<
            DeterministicNeighborDispersalKernel,
            std::default_random_engine>(
            DeterministicNeighborDispersalKernel(Direction::E));
    DynamicDispersalKernel<std::default_random_engine> kernel3(
        kernel2, kernel2, false, 1);

    if (argc != 5) {
        std::cerr
            << "Usage: " << argv[0]
            << " read|stats|write|trips|trace CONFIG_FILE NODE_FILE SEGMENT_FILE\n";
        return 1;
    }
    std::string command = argv[1];
    bool show_stats = false;
    bool write_network = false;
    bool trips = false;
    bool trace = false;
    if (command == "dynamic") {
        show_stats = true;
    }
    else {
        std::cerr << "Unknown sub-command: " << command << "\n";
        std::cerr << "Supported sub-commands are: read, stats, write, trips, trace\n";
        return 1;
    }

    std::string config_file{argv[2]};
    std::ifstream config_stream{config_file};
    if (!config_stream.is_open()) {
        std::cerr << "Failed to open config file: " << config_file << "\n";
        return 1;
    }
    auto raw = read_config(config_stream);
    Config config;

    // config.natural_direction = raw.get("natural_direction",
    // config.natural_direction);
    get_or_use(raw, config.natural_direction, "natural_direction");
    get_or_use(raw, config.natural_kappa, "natural_kappa");
    get_or_use(raw, config.natural_kernel_type, "natural_kernel_type");
    get_or_use(raw, config.natural_scale, "natural_scale");
    get_or_use(raw, config.percent_natural_dispersal, "percent_natural_dispersal");

    dispersers = reproductive_rate * infected;
    int ret = 0;
    for (int i = 0; i < 10000; ++i) {

        auto kernel4 = create_dynamic_natural_kernel<std::default_random_engine>(
            config, dispersers, Network<int>::null_network());
        DynamicDispersalKernel<std::default_random_engine> kernel_d(
            kernel4, nullptr, false, 1);

        simulation.disperse_and_infect(
            i,
            dispersers,
            susceptible,
            exposed,
            infected,
            mortality_tracker,
            total_hosts,
            total_exposed,
            outside_dispersers,
            weather,
            weather_coefficient,
            kernel,
            suitable_cells);
    }

    return ret;
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
    config.use_quarantine = false;
    config.use_spreadrates = true;
    config.spreadrate_frequency = "year";
    config.spreadrate_frequency_n = 1;

    config.set_date_start(2020, 1, 1);
    config.set_date_end(2021, 12, 31);
    config.set_step_unit(StepUnit::Month);
    config.set_step_num_units(1);
    config.use_mortality = false;
    config.mortality_frequency = "year";
    config.mortality_frequency_n = 1;

    config.create_schedules();

    config.deterministic = false;

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
    DoubleRaster weather_coefficient{rows, cols, 0};
    std::vector<std::vector<int>> suitable_cells =
        find_suitable_cells<int>(susceptible);

    IntRaster dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    std::vector<IntRaster> exposed(
        config.latency_period_steps + 1,
        IntRaster(infected.rows(), infected.cols(), 0));

    Simulation<IntRaster, DoubleRaster, int, std::default_random_engine> simulation(
        config.random_seed,
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
            susceptible,
            exposed,
            infected,
            mortality_tracker,
            total_hosts,
            total_exposed,
            outside_dispersers,
            config.weather,
            weather_coefficient,
            kernel,
            suitable_cells);
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
    config.use_quarantine = false;
    config.use_spreadrates = true;
    config.spreadrate_frequency = "year";
    config.spreadrate_frequency_n = 1;

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

    config.deterministic = false;

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
    unsigned rate_num_steps =
        get_number_of_scheduled_actions(config.spread_rate_schedule());
    SpreadRate<Raster<int>> spread_rate(
        infected, config.ew_res, config.ns_res, rate_num_steps, suitable_cells);
    QuarantineEscape<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, 0, suitable_cells);

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
            total_exposed,
            empty_integer,
            mortality_tracker,
            died,
            empty_float,
            empty_float[0],
            treatments,
            zeros,
            outside_dispersers,
            spread_rate,
            quarantine,
            zeros,
            movements,
            Network<int>::null_network(),
            suitable_cells);
    }
    return 0;
}

using SimpleDDKernel = NaturalAnthropogenicDispersalKernel<
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
using DoubleRadialKernel = NaturalAnthropogenicDispersalKernel<
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

int main(int argc, char** argv)
{
    int ret = 0;

    ret += kernel_type_test(argc, argv);

    return ret;
}

#endif  // POPS_TEST
