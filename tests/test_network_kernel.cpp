#include <fstream>
#include <regex>
#include <random>

#include <pops/model.hpp>

using namespace pops;

#include <numeric>

template<typename Number, typename Container>
Number sum(const Container& container)
{
    return std::accumulate(container.begin(), container.end(), Number(0));
}

template<typename Number>
Number sum(const pops::Raster<Number>& raster)
{
    Number ret{0};
    raster.for_each([&ret](const Number& val) { ret += val; });
    return ret;
}

/**
 *
 * Cell center coordinates for bbox 0,100 in both directions and 3x3 raster:
 *
 * ```
 * 100   16.7 50.0 83.3 (x)
 * ----  ---- ---- ----
 * 83.3 | __ | __ | __ |
 * 50.0 | __ | __ | __ |
 * 16.7 |    |    |    |
 * ---   ---- ---- ----
 * 0   0                100
 * (y)
 * ```
 *
 * Nodes (n) and segments (s):
 *
 * ```
 *       16.7 50.0 83.3
 * ----  ---- ---- ----
 * 83.3 | __ | n_ | n_ |  0
 * 50.0 | __ | s_ | n_ |  1
 * 16.7 | n  | s  |    |  2
 * ---   ---- ---- ----
 *        0    1    2
 * ```
 *
 */
int test_model_with_network()
{
    // Data
    Raster<int> infected = {{0, 50, 0}, {0, 0, 50}, {0, 0, 0}};
    Raster<int> susceptible = {{100, 100, 100}, {100, 0, 100}, {100, 0, 100}};
    auto total_hosts = infected + susceptible;
    Raster<int> total_populations = {{100, 100, 100}, {100, 100, 100}, {100, 100, 100}};
    Raster<int> total_exposed(infected.rows(), infected.cols(), 0);
    // Reference data (to be modified later)
    auto original_infected = infected;
    auto expected_susceptible = susceptible;
    // Simulation data
    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;
    // Empty data
    Raster<int> zeros(infected.rows(), infected.cols(), 0);
    std::vector<Raster<double>> empty_floats;
    std::vector<Raster<int>> empty_ints(
        1, Raster<int>(infected.rows(), infected.cols(), 0));
    // Config
    Config config;
    config.random_seed = 0;
    config.model_type = "SI";
    config.reproductive_rate = 100;
    config.natural_kernel_type = "cauchy";
    config.natural_scale = 0.1;
    config.anthro_kernel_type = "network";
    config.network_min_time = 1;
    config.network_max_time = 9;
    config.network_speed = 33.3;
    config.use_anthropogenic_kernel = true;
    config.percent_natural_dispersal = 0;
    config.anthro_scale = config.natural_scale;  // Unused, but we need to set it.
    config.rows = 3;
    config.cols = 3;
    config.ew_res = 33.3;
    config.ns_res = 33.3;

    config.set_date_start(2001, 3, 1);
    config.set_date_end(2001, 3, 3);
    config.create_schedules();

    config.bbox.north = 100;
    config.bbox.south = 0;
    config.bbox.east = 100;
    config.bbox.west = 0;

    Network<int> network{
        config.bbox, config.ew_res, config.ns_res, config.network_speed};
    std::stringstream node_stream{
        "1,16.7,16.7\n"
        "2,50.0,83.3\n"
        "3,83.3,83.3\n"
        "4,83.3,50.0\n"};
    std::stringstream segment_stream{
        "1,2,16.7;16.7;50.0;16.7;50.0;50.0;50.0;83.3\n"
        "4,3,83.3;50.0;83.3;83.3\n"};
    network.load(node_stream, segment_stream);

    // Objects
    std::vector<std::vector<int>> suitable_cells =
        find_suitable_cells<int>(total_hosts);
    Treatments<Raster<int>, Raster<double>> treatments(config.scheduler());
    SpreadRate<Raster<int>> spread_rate(
        infected, config.ew_res, config.ns_res, 0, suitable_cells);
    QuarantineEscape<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, 0, suitable_cells);
    std::vector<std::vector<int>> movements;
    Model<Raster<int>, Raster<double>, Raster<double>::IndexType> model(config);
    // Run
    for (unsigned step = 0; step < config.scheduler().get_num_steps(); ++step) {
        model.run_step(
            step,
            infected,
            susceptible,
            total_populations,
            total_hosts,
            dispersers,
            total_exposed,
            empty_ints,
            empty_ints,
            zeros,
            empty_floats,
            empty_floats[0],
            treatments,
            zeros,
            outside_dispersers,
            spread_rate,
            quarantine,
            zeros,
            movements,
            network,
            suitable_cells);
    }

    int ret = 0;
    std::vector<std::pair<int, int>> should_be_same{{0, 0}, {1, 0}, {2, 2}};
    for (const auto& coords : should_be_same) {
        if (original_infected(coords.first, coords.second)
            != infected(coords.first, coords.second)) {
            std::cerr << "Infected at: " << coords.first << ", " << coords.second
                      << " is different but should be the same"
                      << " (is " << original_infected(coords.first, coords.second)
                      << " but should be " << infected(coords.first, coords.second)
                      << ").\n";
            ret += 1;
        }
    }
    if (sum(original_infected) >= sum(infected)) {
        std::cerr << "New infected not higher than original.\n";
        ret += 1;
    }
    if (ret) {
        std::cerr << "Unexpected (new) infected: \n" << infected << "\n";
        std::cerr << "Original (starting) infected: \n" << original_infected << "\n";
    }
    return ret;
}

int run_tests()
{
    int ret = 0;

    ret += test_model_with_network();

    if (ret)
        std::cerr << "Number of errors in the network kernel test: " << ret << "\n";
    return ret;
}

int main()
{
    return run_tests();
}
