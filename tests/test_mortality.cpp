#ifdef POPS_TEST

/*
 * Simple compilation test for the PoPS Simulation class.
 *
 * Copyright (C) 2018 by the authors.
 *
 * Authors: Vaclav Petras <wenzeslaus gmail com>
 *          Chris Jones <cmjone25 gmail com>
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

int test_mortality()
{
    Raster<int> infected = {{5, 0}, {0, 0}};
    Raster<int> total_hosts = {{10, 5}, {5, 3}};
    std::vector<Raster<int>> mortality_tracker = {{{3, 0}, {0, 0}}, {{2, 0}, {0, 0}}};
    Raster<int> died = {{0, 0}, {0, 0}};
    Raster<int> expected_died = {{4, 0}, {0, 0}};
    Raster<int> expected_infected = {{1, 0}, {0, 0}};
    Raster<int> expected_total_hosts = {{6, 5}, {5, 3}};
    double mortality_rate = 0.50;
    int mortality_time_lag = 0;
    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    Simulation<Raster<int>, Raster<double>> simulation(
        42, infected.rows(), infected.cols());
    simulation.mortality(
        infected,
        total_hosts,
        mortality_rate,
        mortality_time_lag,
        died,
        mortality_tracker,
        suitable_cells);
    if (died != expected_died) {
        cout << "died (actual, expected):\n"
             << died << "  !=\n"
             << expected_died << "\n";
        return 1;
    }
    if (infected != expected_infected) {
        cout << "infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    if (total_hosts != expected_total_hosts) {
        cout << "total_hosts (actual, expected):\n"
             << total_hosts << "  !=\n"
             << expected_total_hosts << "\n";
        return 1;
    }
    return 0;
}

int main()
{
    int num_errors = 0;

    num_errors += test_mortality();
    std::cout << "Mortality number of errors: " << num_errors << std::endl;
    return num_errors;
}
#endif  // POPS_TEST
