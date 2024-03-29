#ifdef POPS_TEST

/*
 * Simple compilation test for the PoPS Mortality module.
 *
 * Copyright (C) 2018 by the authors.
 *
 * Authors: Chris Jones <cmjone25 gmail com>
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
#include <pops/generator_provider.hpp>

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

int test_move_all_no_exposed()
{
    Raster<int> infected = {{5, 0}, {0, 0}};
    Raster<int> susceptible = {{10, 5}, {5, 3}};
    Raster<int> total_exposed = {{0, 0}, {0, 0}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    auto total_hosts = infected + susceptible + total_exposed + resistant;
    std::vector<Raster<int>> exposed = {{{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}};
    std::vector<Raster<int>> mortality_tracker = {{{3, 0}, {0, 0}}, {{2, 0}, {0, 0}}};
    unsigned step = 1;
    int last_index = 0;
    std::vector<std::vector<int>> movements = {{0, 0, 1, 1, 15}, {0, 1, 1, 0, 5}};
    std::vector<unsigned> movement_schedule = {1, 1};
    // set expected new values after movement
    Raster<int> expected_infected = {{0, 0}, {0, 5}};
    Raster<int> expected_susceptible = {{0, 0}, {10, 13}};
    Raster<int> expected_total_exposed = {{0, 0}, {0, 0}};
    Raster<int> expected_resistant = {{0, 0}, {0, 0}};
    Raster<int> expected_total_hosts = {{0, 0}, {10, 18}};
    std::vector<Raster<int>> expected_exposed = {{{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}};
    std::vector<Raster<int>> expected_mortality_tracker = {
        {{0, 0}, {0, 3}}, {{0, 0}, {0, 2}}};
    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    DefaultSingleGeneratorProvider generator(42);
    Simulation<Raster<int>, Raster<double>> simulation(
        infected.rows(), infected.cols());
    simulation.movement(
        infected,
        susceptible,
        mortality_tracker,
        exposed,
        resistant,
        total_hosts,
        total_exposed,
        step,
        last_index,
        movements,
        movement_schedule,
        suitable_cells,
        generator);
    int errors{0};
    if (infected != expected_infected) {
        cout << "test_move_all_no_exposed: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        ++errors;
    }
    if (susceptible != expected_susceptible) {
        cout << "test_move_all_no_exposed: susceptible (actual, expected):\n"
             << susceptible << "  !=\n"
             << expected_susceptible << "\n";
        ++errors;
    }
    if (total_exposed != expected_total_exposed) {
        cout << "test_move_all_no_exposed: total_exposed (actual, expected):\n"
             << total_exposed << "  !=\n"
             << expected_total_exposed << "\n";
        ++errors;
    }
    if (resistant != expected_resistant) {
        cout << "test_move_all_no_exposed: resistant (actual, expected):\n"
             << resistant << "  !=\n"
             << expected_resistant << "\n";
        ++errors;
    }
    if (total_hosts != expected_total_hosts) {
        cout << "test_move_all_no_exposed: total_hosts (actual, expected):\n"
             << total_hosts << "  !=\n"
             << expected_total_hosts << "\n";
        ++errors;
    }
    if (exposed[0] != expected_exposed[0]) {
        cout << "test_move_all_no_exposed: exposeds (actual, expected):\n"
             << exposed[0] << "  !=\n"
             << expected_exposed[0] << "\n";
        ++errors;
    }
    if (mortality_tracker[0] != expected_mortality_tracker[0]) {
        cout << "test_move_all_no_exposed: mortality_tracker (actual, expected):\n"
             << mortality_tracker[0] << "  !=\n"
             << expected_mortality_tracker[0] << "\n";
        ++errors;
    }
    if (exposed[1] != expected_exposed[1]) {
        cout << "test_move_all_no_exposed: exposeds (actual, expected):\n"
             << exposed[1] << "  !=\n"
             << expected_exposed[1] << "\n";
        ++errors;
    }
    if (mortality_tracker[1] != expected_mortality_tracker[1]) {
        cout << "test_move_all_no_exposed: mortality_tracker (actual, expected):\n"
             << mortality_tracker[1] << "  !=\n"
             << expected_mortality_tracker[1] << "\n";
        ++errors;
    }
    return errors;
}

int test_move_all_exposed()
{
    Raster<int> infected = {{5, 0}, {0, 0}};
    Raster<int> susceptible = {{10, 5}, {5, 3}};
    Raster<int> total_exposed = {{5, 0}, {0, 0}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    auto total_hosts = infected + susceptible + total_exposed + resistant;
    std::vector<Raster<int>> exposed = {{{2, 0}, {0, 0}}, {{3, 0}, {0, 0}}};
    std::vector<Raster<int>> mortality_tracker = {{{3, 0}, {0, 0}}, {{2, 0}, {0, 0}}};
    unsigned step = 1;
    int last_index = 0;
    std::vector<std::vector<int>> movements = {{0, 0, 1, 1, 20}, {0, 1, 1, 0, 5}};
    std::vector<unsigned> movement_schedule = {1, 1};
    // set expected new values after movement
    Raster<int> expected_infected = {{0, 0}, {0, 5}};
    Raster<int> expected_susceptible = {{0, 0}, {10, 13}};
    Raster<int> expected_total_exposed = {{0, 0}, {0, 5}};
    Raster<int> expected_resistant = {{0, 0}, {0, 0}};
    Raster<int> expected_total_hosts = {{0, 0}, {10, 23}};
    std::vector<Raster<int>> expected_exposed = {{{0, 0}, {0, 2}}, {{0, 0}, {0, 3}}};
    std::vector<Raster<int>> expected_mortality_tracker = {
        {{0, 0}, {0, 3}}, {{0, 0}, {0, 2}}};
    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    DefaultSingleGeneratorProvider generator(42);
    Simulation<Raster<int>, Raster<double>> simulation(
        infected.rows(), infected.cols());
    simulation.movement(
        infected,
        susceptible,
        mortality_tracker,
        exposed,
        resistant,
        total_hosts,
        total_exposed,
        step,
        last_index,
        movements,
        movement_schedule,
        suitable_cells,
        generator);
    int errors{0};
    if (infected != expected_infected) {
        cout << "test_move_all_exposed: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        ++errors;
    }
    if (susceptible != expected_susceptible) {
        cout << "test_move_all_exposed: susceptible (actual, expected):\n"
             << susceptible << "  !=\n"
             << expected_susceptible << "\n";
        ++errors;
    }
    if (total_exposed != expected_total_exposed) {
        cout << "test_move_all_exposed: total_exposed (actual, expected):\n"
             << total_exposed << "  !=\n"
             << expected_total_exposed << "\n";
        ++errors;
    }
    if (resistant != expected_resistant) {
        cout << "test_move_all_exposed: resistant (actual, expected):\n"
             << resistant << "  !=\n"
             << expected_resistant << "\n";
        ++errors;
    }
    if (total_hosts != expected_total_hosts) {
        cout << "test_move_all_exposed: total_hosts (actual, expected):\n"
             << total_hosts << "  !=\n"
             << expected_total_hosts << "\n";
        ++errors;
    }
    if (exposed[0] != expected_exposed[0]) {
        cout << "test_move_all_exposed: exposeds (actual, expected):\n"
             << exposed[0] << "  !=\n"
             << expected_exposed[0] << "\n";
        ++errors;
    }
    if (mortality_tracker[0] != expected_mortality_tracker[0]) {
        cout << "test_move_all_exposed: mortality_tracker (actual, expected):\n"
             << mortality_tracker[0] << "  !=\n"
             << expected_mortality_tracker[0] << "\n";
        ++errors;
    }
    if (exposed[1] != expected_exposed[1]) {
        cout << "test_move_all_exposed: exposeds (actual, expected):\n"
             << exposed[1] << "  !=\n"
             << expected_exposed[1] << "\n";
        ++errors;
    }
    if (mortality_tracker[1] != expected_mortality_tracker[1]) {
        cout << "test_move_all_exposed: mortality_tracker (actual, expected):\n"
             << mortality_tracker[1] << "  !=\n"
             << expected_mortality_tracker[1] << "\n";
        ++errors;
    }
    return errors;
}

int test_move_add_suitable_cell()
{
    Raster<int> infected = {{5, 0}, {0, 0}};
    Raster<int> susceptible = {{10, 5}, {5, 0}};
    Raster<int> total_exposed = {{5, 0}, {0, 0}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    auto total_hosts = infected + susceptible + total_exposed + resistant;
    std::vector<Raster<int>> exposed = {{{2, 0}, {0, 0}}, {{3, 0}, {0, 0}}};
    std::vector<Raster<int>> mortality_tracker = {{{3, 0}, {0, 0}}, {{2, 0}, {0, 0}}};
    unsigned step = 1;
    int last_index = 0;
    std::vector<std::vector<int>> movements = {{0, 0, 1, 1, 20}, {0, 1, 1, 0, 5}};
    std::vector<unsigned> movement_schedule = {1, 1};
    // set expected new values after movement
    Raster<int> expected_infected = {{0, 0}, {0, 5}};
    Raster<int> expected_susceptible = {{0, 0}, {10, 10}};
    Raster<int> expected_total_exposed = {{0, 0}, {0, 5}};
    Raster<int> expected_resistant = {{0, 0}, {0, 0}};
    Raster<int> expected_total_hosts = {{0, 0}, {10, 20}};
    std::vector<Raster<int>> expected_exposed = {{{0, 0}, {0, 2}}, {{0, 0}, {0, 3}}};
    std::vector<Raster<int>> expected_mortality_tracker = {
        {{0, 0}, {0, 3}}, {{0, 0}, {0, 2}}};
    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}};
    std::vector<std::vector<int>> expected_suitable_cells = {
        {0, 0}, {0, 1}, {1, 0}, {1, 1}};
    DefaultSingleGeneratorProvider generator(42);
    Simulation<Raster<int>, Raster<double>> simulation(
        infected.rows(), infected.cols());
    simulation.movement(
        infected,
        susceptible,
        mortality_tracker,
        exposed,
        resistant,
        total_hosts,
        total_exposed,
        step,
        last_index,
        movements,
        movement_schedule,
        suitable_cells,
        generator);
    int errors{0};
    if (infected != expected_infected) {
        cout << "test_move_add_suitable_cell: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        ++errors;
    }
    if (susceptible != expected_susceptible) {
        cout << "test_move_add_suitable_cell: susceptible (actual, expected):\n"
             << susceptible << "  !=\n"
             << expected_susceptible << "\n";
        ++errors;
    }
    if (total_exposed != expected_total_exposed) {
        cout << "test_move_add_suitable_cell: total_exposed (actual, expected):\n"
             << total_exposed << "  !=\n"
             << expected_total_exposed << "\n";
        ++errors;
    }
    if (resistant != expected_resistant) {
        cout << "test_move_add_suitable_cell: resistant (actual, expected):\n"
             << resistant << "  !=\n"
             << expected_resistant << "\n";
        ++errors;
    }
    if (total_hosts != expected_total_hosts) {
        cout << "test_move_add_suitable_cell: total_hosts (actual, expected):\n"
             << total_hosts << "  !=\n"
             << expected_total_hosts << "\n";
        ++errors;
    }
    if (exposed[0] != expected_exposed[0]) {
        cout << "test_move_add_suitable_cell: exposeds (actual, expected):\n"
             << exposed[0] << "  !=\n"
             << expected_exposed[0] << "\n";
        ++errors;
    }
    if (mortality_tracker[0] != expected_mortality_tracker[0]) {
        cout << "test_move_add_suitable_cell: mortality_tracker (actual, expected):\n"
             << mortality_tracker[0] << "  !=\n"
             << expected_mortality_tracker[0] << "\n";
        ++errors;
    }
    if (exposed[1] != expected_exposed[1]) {
        cout << "test_move_add_suitable_cell: exposeds (actual, expected):\n"
             << exposed[1] << "  !=\n"
             << expected_exposed[1] << "\n";
        ++errors;
    }
    if (mortality_tracker[1] != expected_mortality_tracker[1]) {
        cout << "test_move_add_suitable_cell: mortality_tracker (actual, expected):\n"
             << mortality_tracker[1] << "  !=\n"
             << expected_mortality_tracker[1] << "\n";
        ++errors;
    }
    if (suitable_cells[0][0] != expected_suitable_cells[0][0]) {
        cout << "test_move_add_suitable_cell: mortality_tracker (actual, expected):\n"
             << suitable_cells[0][0] << "  !=\n"
             << expected_suitable_cells[0][0] << "\n";
        ++errors;
    }
    if (suitable_cells[0][1] != expected_suitable_cells[0][1]) {
        cout << "test_move_add_suitable_cell: mortality_tracker (actual, expected):\n"
             << suitable_cells[0][1] << "  !=\n"
             << expected_suitable_cells[0][1] << "\n";
        ++errors;
    }
    if (suitable_cells[1][0] != expected_suitable_cells[1][0]) {
        cout << "test_move_add_suitable_cell: mortality_tracker (actual, expected):\n"
             << suitable_cells[1][0] << "  !=\n"
             << expected_suitable_cells[1][0] << "\n";
        ++errors;
    }
    if (suitable_cells[1][1] != expected_suitable_cells[1][1]) {
        cout << "test_move_add_suitable_cell: mortality_tracker (actual, expected):\n"
             << suitable_cells[1][1] << "  !=\n"
             << expected_suitable_cells[1][1] << "\n";
        ++errors;
    }
    if (suitable_cells[2][0] != expected_suitable_cells[2][0]) {
        cout << "test_move_add_suitable_cell: mortality_tracker (actual, expected):\n"
             << suitable_cells[2][0] << "  !=\n"
             << expected_suitable_cells[2][0] << "\n";
        ++errors;
    }
    if (suitable_cells[2][1] != expected_suitable_cells[2][1]) {
        cout << "test_move_add_suitable_cell: mortality_tracker (actual, expected):\n"
             << suitable_cells[2][1] << "  !=\n"
             << expected_suitable_cells[2][1] << "\n";
        ++errors;
    }
    if (suitable_cells[3][0] != expected_suitable_cells[3][0]) {
        cout << "test_move_add_suitable_cell: mortality_tracker (actual, expected):\n"
             << suitable_cells[3][0] << "  !=\n"
             << expected_suitable_cells[3][0] << "\n";
        ++errors;
    }
    if (suitable_cells[3][1] != expected_suitable_cells[3][1]) {
        cout << "test_move_add_suitable_cell: mortality_tracker (actual, expected):\n"
             << suitable_cells[3][1] << "  !=\n"
             << expected_suitable_cells[3][1] << "\n";
        ++errors;
    }
    return errors;
}

int main()
{
    int num_errors = 0;
    num_errors += test_move_all_no_exposed();
    num_errors += test_move_all_exposed();
    std::cout << "Movement number of errors: " << num_errors << std::endl;
    return num_errors;
}
#endif  // POPS_TEST
