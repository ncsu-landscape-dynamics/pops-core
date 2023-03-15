#ifdef POPS_TEST

/*
 * Simple test for the PoPS survival rate functionality.
 *
 * Copyright (C) 2022 by the authors.
 *
 * Authors: Anna Petrasova <akratoc ncsu edu>
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

#include <pops/raster.hpp>
#include <pops/simulation.hpp>
#include <pops/simple_generator.hpp>

#include <map>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>

using std::string;
using std::cout;
using std::cerr;
using std::endl;

using namespace pops;

int test_survive()
{
    Raster<int> infected = {{5, 0}, {0, 0}};
    Raster<int> susceptible = {{10, 5}, {5, 3}};
    Raster<int> total_exposed = {{10, 0}, {0, 0}};
    std::vector<Raster<int>> exposed = {{{3, 0}, {0, 0}}, {{7, 0}, {0, 0}}};
    std::vector<Raster<int>> mortality_tracker = {{{3, 0}, {0, 0}}, {{2, 0}, {0, 0}}};
    Raster<double> survival_rate = {{0.3, 1}, {0, 0}};
    // set expected new values after movement
    Raster<int> expected_infected = {{2, 0}, {0, 0}};
    Raster<int> expected_susceptible = {{20, 5}, {5, 3}};
    Raster<int> expected_total_exposed = {{3, 0}, {0, 0}};
    // this depends on seed, total should be 3
    std::vector<Raster<int>> expected_exposed = {{{1, 0}, {0, 0}}, {{2, 0}, {0, 0}}};
    // this depends on seed, total should be 2
    std::vector<Raster<int>> expected_mortality_tracker = {
        {{2, 0}, {0, 0}}, {{0, 0}, {0, 0}}};
    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    SimpleGeneratorProvider generator(42);
    Simulation<Raster<int>, Raster<double>> simulation(
        42, infected.rows(), infected.cols());
    simulation.remove_percentage(
        infected,
        susceptible,
        mortality_tracker,
        exposed,
        total_exposed,
        survival_rate,
        suitable_cells,
        generator);
    if (infected != expected_infected) {
        cout << "infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    if (susceptible != expected_susceptible) {
        cout << "susceptible (actual, expected):\n"
             << susceptible << "  !=\n"
             << expected_susceptible << "\n";
        return 1;
    }
    if (total_exposed != expected_total_exposed) {
        cout << "total_exposed (actual, expected):\n"
             << total_exposed << "  !=\n"
             << expected_total_exposed << "\n";
        return 1;
    }

    if (exposed[0] != expected_exposed[0]) {
        cout << "exposeds (actual, expected):\n"
             << exposed[0] << "  !=\n"
             << expected_exposed[0] << "\n";
        return 1;
    }
    if (mortality_tracker[0] != expected_mortality_tracker[0]) {
        cout << "mortality_tracker (actual, expected):\n"
             << mortality_tracker[0] << "  !=\n"
             << expected_mortality_tracker[0] << "\n";
        return 1;
    }
    if (exposed[1] != expected_exposed[1]) {
        cout << "exposeds (actual, expected):\n"
             << exposed[1] << "  !=\n"
             << expected_exposed[1] << "\n";
        return 1;
    }
    if (mortality_tracker[1] != expected_mortality_tracker[1]) {
        cout << "mortality_tracker (actual, expected):\n"
             << mortality_tracker[1] << "  !=\n"
             << expected_mortality_tracker[1] << "\n";
        return 1;
    }
    return 0;
}

int main()
{
    int num_errors = 0;
    num_errors += test_survive();
    std::cout << "Survival rate number of errors: " << num_errors << std::endl;
    return num_errors;
}
#endif  // POPS_TEST
