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
#include <pops/raster.hpp>
#include <pops/neighbor_kernel.hpp>
#include <pops/simulation.hpp>

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

int main()
{
    int num_errors = 0;

    num_errors += test_infected_arrive();
    return num_errors;
}

#endif  // POPS_TEST
