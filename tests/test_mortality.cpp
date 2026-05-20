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
        infected.rows(), infected.cols());
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

void warn_about_tracker_size(
    const std::vector<Raster<int>>& mortality_tracker,
    double mortality_rate,
    int mortality_time_lag)
{
    if (mortality_tracker.size() != 1 / mortality_rate + mortality_time_lag) {
        cerr << "Non-enforced requirement of mortality tracker vector size "
                "broken in the mortality time lag test: "
             << "actual: " << mortality_tracker.size()
             << ", expected: " << 1 / mortality_rate + mortality_time_lag << "\n";
    }
}

int compare_lists_equal(
    const std::vector<Raster<int>>& reference,
    const std::vector<Raster<int>>& actual,
    const std::string& context,
    int step)
{
    int ret = 0;
    if (reference.size() != actual.size()) {
        cerr << context << " (step " << step
             << "): compare_lists_equal: reference and actual differs in size: "
             << reference.size() << " != " << actual.size() << "\n";
        ++ret;
    }
    for (size_t i = 0; i < std::min(reference.size(), actual.size()); ++i) {
        if (reference.at(i) != actual.at(i)) {
            cerr << context << " (step " << step
                 << "): compare_lists_equal: reference and actual differs at index "
                 << i << " (max index is " << reference.size() - 1 << "): \n"
                 << reference.at(i) << " !=\n"
                 << actual.at(i) << "\n";
            ++ret;
        }
    }
    return ret;
}

/**
 * Test mortality with time lag
 */
int test_mortality_lag()
{
    int ret = 0;

    Raster<int> infected = {{6, 0}, {0, 0}};
    Raster<int> total_hosts = {{14, 80}, {6, 4}};
    Raster<int> died = {{0, 0}, {0, 0}};
    Raster<int> expected_died = {{0, 0}, {0, 0}};
    Raster<int> expected_infected = infected;
    Raster<int> expected_total_hosts = total_hosts;
    double mortality_rate = 0.50;
    int mortality_time_lag = 2;
    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};

    Simulation<Raster<int>, Raster<double>> simulation(
        infected.rows(), infected.cols());

    int step = 1;

    Raster<int> newly_infected = total_hosts - infected;
    infected += newly_infected;
    expected_infected += newly_infected;
    std::vector<Raster<int>> mortality_tracker = {
        {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, newly_infected};
    warn_about_tracker_size(mortality_tracker, mortality_rate, mortality_time_lag);

    simulation.mortality(
        infected,
        total_hosts,
        mortality_rate,
        mortality_time_lag,
        died,
        mortality_tracker,
        suitable_cells);
    ret += compare_lists_equal(
        {{{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, newly_infected, {{0, 0}, {0, 0}}},
        mortality_tracker,
        "time lag",
        step);
    if (died != expected_died) {
        cout << "time lag: died (actual, expected), step " << step << ":\n"
             << died << "  !=\n"
             << expected_died << "\n";
        ++ret;
    }
    if (infected != expected_infected) {
        cout << "time lag: infected (actual, expected), step " << step << ":\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        ++ret;
    }
    if (total_hosts != expected_total_hosts) {
        cout << "time lag: total_hosts (actual, expected), step " << step << ":\n"
             << total_hosts << "  !=\n"
             << expected_total_hosts << "\n";
        ++ret;
    }

    ++step;
    simulation.mortality(
        infected,
        total_hosts,
        mortality_rate,
        mortality_time_lag,
        died,
        mortality_tracker,
        suitable_cells);
    ret += compare_lists_equal(
        {{{0, 0}, {0, 0}}, newly_infected, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}},
        mortality_tracker,
        "time lag",
        step);
    if (died != expected_died) {
        cout << "time lag: died (actual, expected), step " << step << ":\n"
             << died << "  !=\n"
             << expected_died << "\n";
        ++ret;
    }
    if (infected != expected_infected) {
        cout << "time lag: infected (actual, expected), step " << step << ":\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        ++ret;
    }
    if (total_hosts != expected_total_hosts) {
        cout << "time lag: total_hosts (actual, expected), step " << step << ":\n"
             << total_hosts << "  !=\n"
             << expected_total_hosts << "\n";
        ++ret;
    }

    ++step;

    expected_died = {{8 / 2, 80 / 2}, {6 / 2, 4 / 2}};  // 6 infected are initial state
    expected_infected -= expected_died;
    expected_total_hosts -= expected_died;
    warn_about_tracker_size(mortality_tracker, mortality_rate, mortality_time_lag);
    simulation.mortality(
        infected,
        total_hosts,
        mortality_rate,
        mortality_time_lag,
        died,
        mortality_tracker,
        suitable_cells);
    ret += compare_lists_equal(
        {newly_infected - expected_died,
         {{0, 0}, {0, 0}},
         {{0, 0}, {0, 0}},
         {{0, 0}, {0, 0}}},
        mortality_tracker,
        "time lag",
        step);
    if (died != expected_died) {
        cout << "time lag: died (actual, expected), step " << step << ":\n"
             << died << "  !=\n"
             << expected_died << "\n";
        ++ret;
    }
    if (infected != expected_infected) {
        cout << "time lag: infected (actual, expected), step " << step << ":\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        ++ret;
    }
    if (total_hosts != expected_total_hosts) {
        cout << "time lag: total_hosts (actual, expected), step " << step << ":\n"
             << total_hosts << "  !=\n"
             << expected_total_hosts << "\n";
        ++ret;
    }

    // We don't check cumulative died, but died in the step.
    died.zero();

    ++step;

    // Second half of the infected.
    expected_died = {{8 / 2, 80 / 2}, {6 / 2, 4 / 2}};
    expected_infected -= expected_died;
    expected_total_hosts -= expected_died;
    warn_about_tracker_size(mortality_tracker, mortality_rate, mortality_time_lag);
    simulation.mortality(
        infected,
        total_hosts,
        mortality_rate,
        mortality_time_lag,
        died,
        mortality_tracker,
        suitable_cells);
    ret += compare_lists_equal(
        {{{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}}},
        mortality_tracker,
        "time lag",
        step);
    if (died != expected_died) {
        cout << "time lag: died (actual, expected), step " << step << ":\n"
             << died << "  !=\n"
             << expected_died << "\n";
        ++ret;
    }
    if (infected != expected_infected) {
        cout << "time lag: infected (actual, expected), step " << step << ":\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        ++ret;
    }
    if (total_hosts != expected_total_hosts) {
        cout << "time lag: total_hosts (actual, expected), step " << step << ":\n"
             << total_hosts << "  !=\n"
             << expected_total_hosts << "\n";
        ++ret;
    }
    return ret;
}

int main()
{
    int num_errors = 0;

    num_errors += test_mortality();
    num_errors += test_mortality_lag();
    std::cout << "Mortality number of errors: " << num_errors << std::endl;
    return num_errors;
}
#endif  // POPS_TEST
