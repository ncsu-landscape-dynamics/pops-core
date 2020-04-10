#ifdef POPS_TEST

/*
 * Simple compilation test for the PoPS Simulation class.
 *
 * Copyright (C) 2018-2020 by the authors.
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

#include "raster.hpp"
#include "neighbor_kernel.hpp"
#include "radial_kernel.hpp"
#include "simulation.hpp"

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

template <typename T>
void print_vector(const std::vector<T>& v)
{
    for (auto i : v) {
        cout << i;
    }
    cout << "\n";
}

int test_rotate_left_by_one(std::vector<int> a, std::vector<int> b)
{
    rotate_left_by_one(a);
    if (a != b) {
        cout << "Rotated vector not correct\n";
        print_vector(a);
        print_vector(b);
        return 1;
    }
    return 0;
}

int test_with_neighbor_kernel()
{
    Raster<int> infected = {{5, 0}, {0, 0}};
    Raster<int> mortality_tracker = {{0, 0}, {0, 0}};
    // Susceptible and total are set in a way that there won't be any
    // dilution effect and the disperser will always establish given the
    // selected random seed. Establishment probability is high and with
    // the given seed we don't get any random numbers in establishment
    // test higher than that. (The weather is disabled.)
    Raster<int> susceptible = {{10, 6}, {14, 15}};
    // add a lot of hosts, so that exposing or infecting them won't
    // chanage the susceptible/total ratio much
    susceptible += 100000;
    // we want to minimize the dilution effect
    Raster<int> total_plants = susceptible;
    Raster<double> temperature = {{5, 0}, {0, 0}};
    Raster<double> weather_coefficient = {{0, 0}, {0, 0}};

    Raster<int> expected_exposed = {{0, 10}, {0, 0}};
    Raster<int> expected_mortality_tracker = expected_exposed;

    Raster<int> exposed(infected.rows(), infected.cols(), 0);
    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;
    bool weather = false;
    double reproductive_rate = 2;
    unsigned latency_period_steps = 2;
    DeterministicNeighborDispersalKernel kernel(Direction::E);
    Simulation<Raster<int>, Raster<double>> simulation(
                model_type_from_string("SI"),
                latency_period_steps,
                42,
                infected.rows(),
                infected.cols()
                );
    dispersers = reproductive_rate * infected;
    // cout << dispersers;
    simulation.disperse(dispersers, susceptible, exposed,
                        mortality_tracker, total_plants,
                        outside_dispersers, weather, weather_coefficient,
                        kernel);
    if (!outside_dispersers.empty()) {
        cout << "There are outside_dispersers (" << outside_dispersers.size() << ") but there should be none\n";
        return 1;
    }
    if (exposed != expected_exposed) {
        cout << "Neighbor kernel test (actual, expected):\n" << exposed << "  !=\n" << expected_exposed << "\n";
        return 1;
    }
    if (mortality_tracker != expected_mortality_tracker) {
        cout << "Neighbor kernel test mortality tracker (actual, expected):\n" << mortality_tracker << "  !=\n" << expected_mortality_tracker << "\n";
        return 1;
    }
    return 0;
}

int disperse_and_infect_postcondition(int step, const std::vector<Raster<int>>& exposed)
{
    Raster<int> zeros(exposed[0].rows(), exposed[0].cols(), 0);
    if (exposed.back() != zeros) {
        cout << "SEI: disperse_and_infect post-condition not met in step " << step << "\n";
        return 1;
    }
    return 0;
}

int exposed_state(int step, const std::vector<Raster<int>>& exposed, const Raster<int>& expected_exposed)
{
    int ret = 0;
    Raster<int> zeros(exposed[0].rows(), exposed[0].cols(), 0);
    for (int i = 0; i < exposed.size(); ++i) {
        if (i >= int(exposed.size()) - step - 2 && i < exposed.size() - 1) {
            if (exposed[i] != expected_exposed) {
                cout << "SEI test exposed[" << i << "] (actual, expected):\n" << exposed[i] << "  !=\n" << expected_exposed << "\n";
                print_vector(exposed);
                ret += 1;
            }
        }
        else {
            if (exposed[i] != zeros) {
                cout << "SEI test exposed[" << i << "] (actual, expected zeros):\n" << exposed[i] << "\n";
                print_vector(exposed);
                ret += 1;
            }
        }
    }
    return ret;
}

int test_with_sei()
{
    Raster<int> infected = {{5, 0}, {0, 0}};
    Raster<int> mortality_tracker = {{0, 0}, {0, 0}};
    // Susceptible and total are set in a way that there won't be any
    // dilution effect and the disperser will always establish given the
    // selected random seed. Establishment probability is high and with
    // the given seed we don't get any random numbers in establishment
    // test higher than that. (The weather is disabled.)
    Raster<int> susceptible = {{10, 6}, {14, 15}};
    // add a lot of hosts, so that exposing or infecting them won't
    // chanage the susceptible/total ratio much
    susceptible += 100000;
    // we want to minimize the dilution effect
    Raster<int> total_plants = susceptible;
    Raster<double> temperature = {{5, 0}, {0, 0}};
    Raster<double> weather_coefficient = {{0, 0}, {0, 0}};
    Raster<int> zeros(infected.rows(), infected.cols(), 0);

    Raster<int> expected_infected = infected;
    Raster<int> expected_exposed = {{0, 10}, {0, 0}};

    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;
    bool weather = false;
    double reproductive_rate = 2;
    unsigned latency_period_steps = 3;

    std::vector<Raster<int>> exposed(
                latency_period_steps + 1,
                Raster<int>(infected.rows(), infected.cols(), 0));

    DeterministicNeighborDispersalKernel kernel(Direction::E);
    Simulation<Raster<int>, Raster<double>> simulation(
                model_type_from_string("SEI"),
                latency_period_steps,
                42,
                infected.rows(),
                infected.cols()
                );
    dispersers = reproductive_rate * infected;
    int step = 0;
    int ret = 0;
    simulation.disperse_and_infect(
                step, dispersers, susceptible,
                exposed, infected,
                mortality_tracker, total_plants,
                outside_dispersers, weather,
                weather_coefficient,
                kernel);
    if (infected != expected_infected) {
        cout << "SEI test infected (actual, expected):\n" << infected << "  !=\n" << expected_infected << "\n";
        ret += 1;
    }
    if (mortality_tracker != zeros) {
        cout << "SEI test mortality tracker (actual, expected zeros):\n" << mortality_tracker << "\n";
        ret += 1;
    }
    print_vector(exposed);
    cout << infected << "\n\n";
    ret += disperse_and_infect_postcondition(step, exposed);
    simulation.disperse_and_infect(
                ++step, dispersers, susceptible,
                exposed, infected,
                mortality_tracker, total_plants,
                outside_dispersers, weather,
                weather_coefficient,
                kernel);
    print_vector(exposed);
    cout << infected << "\n\n";
    ret += disperse_and_infect_postcondition(step, exposed);
    simulation.disperse_and_infect(
                ++step, dispersers, susceptible,
                exposed, infected,
                mortality_tracker, total_plants,
                outside_dispersers, weather,
                weather_coefficient,
                kernel);
    print_vector(exposed);
    cout << infected << "\n\n";
    ret += disperse_and_infect_postcondition(step, exposed);
    if (!outside_dispersers.empty()) {
        cout << "SEI test: There are outside_dispersers (" << outside_dispersers.size() << ") but there should be none\n";
        ret += 1;
    }
    exposed_state(step, exposed, expected_exposed);
    if (mortality_tracker != zeros) {
        cout << "SEI test mortality tracker (actual, expected zeros):\n" << mortality_tracker << "\n";
        ret += 1;
    }
    simulation.disperse_and_infect(
                ++step, dispersers, susceptible,
                exposed, infected,
                mortality_tracker, total_plants,
                outside_dispersers, weather,
                weather_coefficient,
                kernel);
    print_vector(exposed);
    cout << infected << "\n\n";
    ret += disperse_and_infect_postcondition(step, exposed);
    expected_infected = expected_infected + expected_exposed;
    Raster<int> expected_mortality_tracker = expected_exposed;
    if (!outside_dispersers.empty()) {
        cout << "SEI test: There are outside_dispersers (" << outside_dispersers.size() << ") but there should be none\n";
        ret += 1;
    }
    if (infected != expected_infected) {
        cout << "SEI test infected (actual, expected):\n" << infected << "  !=\n" << expected_infected << "\n";
        ret += 1;
    }
    if (mortality_tracker != expected_mortality_tracker) {
        cout << "SEI test mortality tracker (actual, expected):\n" << mortality_tracker << "  !=\n" << expected_mortality_tracker << "\n";
        ret += 1;
    }
    exposed_state(step, exposed, expected_exposed);
    simulation.disperse_and_infect(
                ++step, dispersers, susceptible,
                exposed, infected,
                mortality_tracker, total_plants,
                outside_dispersers, weather,
                weather_coefficient,
                kernel);
    print_vector(exposed);
    cout << infected << "\n\n";
    ret += disperse_and_infect_postcondition(step, exposed);

    for (int i = 0; i < 10; ++i) {
        simulation.disperse_and_infect(
                    ++step, dispersers, susceptible,
                    exposed, infected,
                    mortality_tracker, total_plants,
                    outside_dispersers, weather,
                    weather_coefficient,
                    kernel);
        print_vector(exposed);
        cout << infected << "\n\n";
        ret += disperse_and_infect_postcondition(step, exposed);
    }

    return ret;
}

int test_calling_all_functions()
{
    Raster<int> infected = {{5, 0}, {0, 0}};
    Raster<int> mortality_tracker = {{0, 0}, {0, 0}};
    Raster<int> susceptible = {{10, 15}, {14, 15}};
    Raster<int> total_plants = {{15, 15}, {14, 15}};
    Raster<double> temperature = {{5, 0}, {0, 0}};
    Raster<double> weather_coefficient = {{0.6, 0.8}, {0.2, 0.8}};
    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;
    DispersalKernelType dispersal_kernel = DispersalKernelType::Cauchy;
    bool weather = true;
    double lethal_temperature = -4.5;
    double reproductive_rate = 4.5;
    double short_distance_scale = 0.0;
    int ew_res = 30;
    int ns_res = 30;
    unsigned step = 1;
    unsigned latency_period_steps = 2;
    unsigned last_index = 0;
    int seed = 42;
    std::vector<std::vector<int>> movements = {{0, 0, 1, 1, 2}, {0, 1, 0, 0, 3}};
    std::vector<unsigned> movement_schedule = {1, 1};
    Simulation<Raster<int>, Raster<double>> simulation(
                model_type_from_string("SI"),
                latency_period_steps,
                seed,
                infected.rows(),
                infected.cols());
    simulation.remove(infected, susceptible, temperature, lethal_temperature);
    simulation.generate(dispersers, infected, weather, weather_coefficient, reproductive_rate);
    RadialDispersalKernel kernel(ew_res, ns_res, dispersal_kernel,
                                 short_distance_scale);
    simulation.movement(infected, susceptible, mortality_tracker, total_plants, step, 
                        last_index, movements, movement_schedule);
    simulation.disperse(dispersers, susceptible, infected,
                        mortality_tracker, total_plants,
                        outside_dispersers, weather, weather_coefficient,
                        kernel);
    cout << "outside_dispersers: " << outside_dispersers.size() << endl;
    return 0;
}

int main()
{
    int ret = 0;

    ret += test_rotate_left_by_one({2}, {2});
    ret += test_rotate_left_by_one({1, 2}, {2, 1});
    ret += test_rotate_left_by_one({1, 2, 3, 4, 5}, {2, 3, 4, 5, 1});

    ret += test_calling_all_functions();
    ret += test_with_neighbor_kernel();
    ret += test_with_sei();

    return ret;
}

#endif  // POPS_TEST
