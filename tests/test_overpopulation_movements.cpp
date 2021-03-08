#ifdef POPS_TEST

/*
 * Test of overpopulation-based movement of pest.
 *
 * Copyright (C) 2021 by the authors.
 *
 * Authors: Vaclav Petras <wenzeslaus gmail com>
 *
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

#include <vector>
#include <pops/model.hpp>
#include <pops/simulation.hpp>
#include <pops/neighbor_kernel.hpp>
#include <pops/raster.hpp>

using namespace pops;

std::vector<std::vector<int>> get_suitable_cells(Raster<int> raster)
{
    std::vector<std::vector<int>> cells;
    for (int row = 0; row < raster.rows(); ++row) {
        for (int col = 0; col < raster.cols(); ++col) {
            if (raster(row, col) > 0) {
                cells.push_back({row, col});
            }
        }
    }
    return cells;
}

int test_infected_arrive()
{
    Raster<int> infected = {{16, 0}, {0, 0}};
    Raster<int> susceptible = {{4, 10}, {20, 15}};
    Raster<int> total_hosts = infected + susceptible;
    std::vector<std::tuple<int, int>> outside_dispersers;
    int seed = 42;
    std::vector<std::vector<int>> suitable_cells = get_suitable_cells(total_hosts);
    Simulation<Raster<int>, Raster<double>> simulation(
        seed, infected.rows(), infected.cols());
    DeterministicNeighborDispersalKernel kernel{Direction::E};
    double overpopulation_percentage = 0.75;
    double leaving_percentage = 0.5;
    simulation.move_overpopulated_pests(
        susceptible,
        infected,
        total_hosts,
        outside_dispersers,
        kernel,
        suitable_cells,
        overpopulation_percentage,
        leaving_percentage);
    int ret = 0;
    Raster<int> expected_infected = {{8, 8}, {0, 0}};
    if (expected_infected != infected) {
        std::cerr << "Unexpected infected: \n" << infected << "\n";
        ret += 1;
    }
    if (outside_dispersers.size() != 0) {
        std::cerr << "Expecting no outside disprersers: " << outside_dispersers.size()
                  << "\n";
        ret += 1;
    }
    simulation.move_overpopulated_pests(
        susceptible,
        infected,
        total_hosts,
        outside_dispersers,
        kernel,
        suitable_cells,
        overpopulation_percentage,
        leaving_percentage);
    expected_infected = {{8, 4}, {0, 0}};
    if (expected_infected != infected) {
        std::cerr << "Unexpected infected: \n" << infected << "\n";
        ret += 1;
    }
    size_t expected_outside_dispersers = 4;
    if (outside_dispersers.size() != expected_outside_dispersers) {
        std::cerr << "Expecting " << expected_outside_dispersers
                  << " outside dispersers, got " << outside_dispersers.size() << "\n";
        ret += 1;
    }
    return ret;
}

int test_model()
{
    // Data
    Raster<int> infected = {{60, 0, 0}, {0, 10, 0}, {0, 0, 2}};
    Raster<int> susceptible = {{40, 100, 100}, {100, 90, 100}, {100, 0, 98}};
    auto total_hosts = infected + susceptible;
    // Reference data (to be modified later)
    auto expected_infected = infected;
    auto expected_susceptible = susceptible;
    // Simulation data
    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;
    // Empty data
    Raster<int> zeros(infected.rows(), infected.cols(), 0);
    std::vector<Raster<double>> empty_floats;
    std::vector<Raster<int>> empty_ints;
    // Config
    Config config;
    config.random_seed = 0;
    config.model_type = "SI";
    config.natural_kernel_type = "cauchy";
    config.natural_scale = 0.1;
    // We are not using anthropo kernel in standard disperser spread, but it is used by
    // overpopulation movements.
    config.use_anthropogenic_kernel = false;
    config.anthro_kernel_type = "cauchy";
    config.anthro_scale = 0.1;
    config.ew_res = 30;
    config.ns_res = 30;
    config.use_overpopulation_movements = true;
    config.overpopulation_percentage = 0.5;
    config.leaving_percentage = 0.75;
    config.create_schedules();
    // More reference data
    auto leaving = infected(0, 0) * config.leaving_percentage;
    // Objects
    const std::vector<std::vector<int>> suitable_cells =
        get_suitable_cells(total_hosts);
    Treatments<Raster<int>, Raster<double>> treatments(config.scheduler());
    SpreadRate<Raster<int>> spread_rate(
        infected, config.ew_res, config.ns_res, 0, suitable_cells);
    QuarantineEscape<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, 0, suitable_cells);
    std::vector<std::vector<int>> movements;
    Model<Raster<int>, Raster<double>, Raster<double>::IndexType> model(config);
    // Run
    model.run_step(
        1,
        infected,
        susceptible,
        total_hosts,
        dispersers,
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
        suitable_cells);
    // Test
    int ret = 0;
    // Test that infected decrease
    expected_infected(0, 0) -= leaving;
    if (expected_infected != infected) {
        std::cerr << "Unexpected infected: \n" << infected << "\n";
        ret += 1;
    }
    // Test that susceptible increase
    expected_susceptible(0, 0) += leaving;
    if (expected_susceptible != susceptible) {
        std::cerr << "Unexpected susceptible: \n" << susceptible << "\n";
        ret += 1;
    }
    // Test that basic counts are correct everywhere
    if (infected + susceptible != total_hosts) {
        std::cerr << "I + S does not match the original total hosts\n";
        std::cerr << "I: \n" << infected << "\n";
        std::cerr << "S: \n" << susceptible << "\n";
        std::cerr << "S + I: \n" << infected + susceptible << "\n";
        std::cerr << "original total: \n" << total_hosts << "\n";
        std::cerr << "diff: \n" << total_hosts - (infected + susceptible) << "\n";
        ret += 1;
    }
    // Test that the leaving ones are in outside
    size_t expected_outside_dispersers = leaving;
    if (outside_dispersers.size() != expected_outside_dispersers) {
        std::cerr << "Expecting " << expected_outside_dispersers
                  << " outside dispersers, got " << outside_dispersers.size() << "\n";
        ret += 1;
    }
    return ret;
}

int main()
{
    int num_errors = 0;

    num_errors += test_infected_arrive();
    num_errors += test_model();
    return num_errors;
}

#endif  // POPS_TEST
