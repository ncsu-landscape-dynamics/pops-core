#ifdef POPS_TEST

/*
 * Simple compilation test for the PoPS deterministic_kernel class.
 *
 * Copyright (C) 2018-2020 by the authors.
 *
 * Authors: Margaret Lawrimore <malawrim ncsu edu>
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

#include <pops/deterministic_kernel.hpp>
#include <pops/cauchy_kernel.hpp>
#include <pops/exponential_kernel.hpp>
#include <pops/exponential_power_kernel.hpp>
#include <pops/lognormal_kernel.hpp>
#include <pops/gamma_kernel.hpp>
#include <pops/simulation.hpp>

using std::string;
using std::cout;

using namespace pops;

int test_with_cauchy_deterministic_kernel()
{
    Raster<int> infected = {{5, 0, 0}, {0, 5, 0}, {0, 0, 2}};
    Raster<int> mortality_tracker = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    Raster<int> susceptible = {{10, 20, 9}, {14, 15, 0}, {3, 0, 2}};
    Raster<int> total_hosts = susceptible;
    Raster<double> temperature = {{5, 0, 0}, {0, 0, 0}, {0, 0, 2}};
    Raster<double> weather_coefficient = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    std::vector<std::vector<int>> movements = {
        {0, 0, 1, 1, 2}, {0, 1, 0, 0, 3}, {0, 1, 1, 0, 2}};
    std::vector<unsigned> movement_schedule = {1, 1};

    Raster<int> expected_mortality_tracker = {{10, 0, 0}, {0, 10, 0}, {0, 0, 2}};
    Raster<int> expected_infected = {{15, 0, 0}, {0, 15, 0}, {0, 0, 4}};

    const std::vector<std::vector<int>> suitable_cell = {
        {0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};

    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;
    bool weather = false;
    double reproductive_rate = 2;
    bool generate_stochasticity = false;
    bool establishment_stochasticity = false;
    bool movement_stochasticity = false;
    // We want everything to establish.
    double establishment_probability = 1;
    // Cauchy
    Simulation<Raster<int>, Raster<double>> simulation(
        42,
        infected.rows(),
        infected.cols(),
        model_type_from_string("SI"),
        0,
        generate_stochasticity,
        establishment_stochasticity,
        movement_stochasticity);
    simulation.generate(
        dispersers,
        infected,
        weather,
        weather_coefficient,
        reproductive_rate,
        suitable_cell);
    auto expected_dispersers = reproductive_rate * infected;
    if (dispersers != expected_dispersers) {
        cout << "Deterministic Kernel Cauchy: dispersers (actual, expected):\n"
             << dispersers << "  !=\n"
             << expected_dispersers << "\n";
        return 1;
    }
    DeterministicDispersalKernel<Raster<int>> deterministicKernel(
        DispersalKernelType::Cauchy, dispersers, 0.9, 30, 30, 0.9);
    // using a smaller scale value since the test raster is so small
    simulation.disperse(
        dispersers,
        susceptible,
        infected,
        mortality_tracker,
        total_hosts,
        outside_dispersers,
        weather,
        weather_coefficient,
        deterministicKernel,
        suitable_cell,
        establishment_probability);
    if (outside_dispersers.size() != 0) {
        cout << "Deterministic Kernel Cauchy: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be 0\n";
        return 1;
    }
    if (infected != expected_infected) {
        cout << "Deterministic Kernel Cauchy: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    if (mortality_tracker != expected_mortality_tracker) {
        cout << "Deterministic Kernel Cauchy: mortality tracker (actual, expected):\n"
             << mortality_tracker << "  !=\n"
             << expected_mortality_tracker << "\n";
        return 1;
    }
    return 0;
}

int test_with_exponential_deterministic_kernel()
{
    Raster<int> infected = {{5, 0, 0}, {0, 5, 0}, {0, 0, 2}};
    Raster<int> mortality_tracker = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    Raster<int> susceptible = {{10, 20, 9}, {14, 15, 0}, {3, 0, 2}};
    Raster<int> total_hosts = susceptible;
    Raster<double> movements = {{0, 0, 1, 1, 2}, {0, 1, 0, 0, 3}, {0, 1, 1, 0, 2}};
    Raster<double> temperature = {{5, 0, 0}, {0, 0, 0}, {0, 0, 2}};
    Raster<double> weather_coefficient = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    std::vector<unsigned> movement_schedule = {1, 1};

    Raster<int> expected_mortality_tracker = {{10, 0, 0}, {0, 10, 0}, {0, 0, 2}};
    Raster<int> expected_infected = {{15, 0, 0}, {0, 15, 0}, {0, 0, 4}};

    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    const std::vector<std::vector<int>> suitable_cell = {
        {0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};

    bool weather = false;
    double reproductive_rate = 2;
    bool generate_stochasticity = false;
    bool establishment_stochasticity = false;
    bool movement_stochasticity = false;
    // We want everything to establish.
    double establishment_probability = 1;
    // Exponential
    Simulation<Raster<int>, Raster<double>> s2(
        42,
        infected.rows(),
        infected.cols(),
        model_type_from_string("SI"),
        0,
        generate_stochasticity,
        establishment_stochasticity,
        movement_stochasticity);
    s2.generate(
        dispersers,
        infected,
        weather,
        weather_coefficient,
        reproductive_rate,
        suitable_cell);
    auto expected_dispersers = reproductive_rate * infected;
    if (dispersers != expected_dispersers) {
        cout << "Deterministic Kernel Exponential: dispersers (actual, expected):\n"
             << dispersers << "  !=\n"
             << expected_dispersers << "\n";
        return 1;
    }
    DeterministicDispersalKernel<Raster<int>> deterministicKernel(
        DispersalKernelType::Exponential, dispersers, 0.99, 30, 30, 1.0);

    s2.disperse(
        dispersers,
        susceptible,
        infected,
        mortality_tracker,
        total_hosts,
        outside_dispersers,
        weather,
        weather_coefficient,
        deterministicKernel,
        suitable_cell,
        establishment_probability);
    if (outside_dispersers.size() != 0) {
        cout << "Deterministic Kernel Exponential: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be 0\n";
        return 1;
    }
    if (infected != expected_infected) {
        cout << "Deterministic Kernel Exponential: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    if (mortality_tracker != expected_mortality_tracker) {
        cout
            << "Deterministic Kernel Exponential: mortality tracker (actual, expected):\n"
            << mortality_tracker << "  !=\n"
            << expected_mortality_tracker << "\n";
        return 1;
    }
    return 0;
}

int test_with_weibull_deterministic_kernel()
{
    Raster<int> infected = {{5, 0, 0}, {0, 5, 0}, {0, 0, 2}};
    Raster<int> mortality_tracker = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    Raster<int> susceptible = {{10, 20, 9}, {14, 15, 0}, {3, 0, 2}};
    Raster<int> total_hosts = susceptible;
    Raster<double> movements = {{0, 0, 1, 1, 2}, {0, 1, 0, 0, 3}, {0, 1, 1, 0, 2}};
    Raster<double> temperature = {{5, 0, 0}, {0, 0, 0}, {0, 0, 2}};
    Raster<double> weather_coefficient = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    std::vector<unsigned> movement_schedule = {1, 1};

    Raster<int> expected_mortality_tracker = {{10, 0, 0}, {0, 10, 0}, {0, 0, 2}};
    Raster<int> expected_infected = {{15, 0, 0}, {0, 15, 0}, {0, 0, 4}};

    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    const std::vector<std::vector<int>> suitable_cell = {
        {0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};

    bool weather = false;
    double reproductive_rate = 2;
    bool generate_stochasticity = false;
    bool establishment_stochasticity = false;
    bool movement_stochasticity = false;
    // We want everything to establish.
    double establishment_probability = 1;
    Simulation<Raster<int>, Raster<double>> s2(
        42,
        infected.rows(),
        infected.cols(),
        model_type_from_string("SI"),
        0,
        generate_stochasticity,
        establishment_stochasticity,
        movement_stochasticity);
    s2.generate(
        dispersers,
        infected,
        weather,
        weather_coefficient,
        reproductive_rate,
        suitable_cell);
    auto expected_dispersers = reproductive_rate * infected;
    if (dispersers != expected_dispersers) {
        cout << "Deterministic Kernel Weibull: dispersers (actual, expected):\n"
             << dispersers << "  !=\n"
             << expected_dispersers << "\n";
        return 1;
    }
    DeterministicDispersalKernel<Raster<int>> deterministicKernel(
        DispersalKernelType::Weibull, dispersers, 0.99, 30, 30, 1.0, 1.0);

    s2.disperse(
        dispersers,
        susceptible,
        infected,
        mortality_tracker,
        total_hosts,
        outside_dispersers,
        weather,
        weather_coefficient,
        deterministicKernel,
        suitable_cell,
        establishment_probability);
    if (outside_dispersers.size() != 0) {
        cout << "Deterministic Kernel Weibull: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be 0\n";
        return 1;
    }
    if (infected != expected_infected) {
        cout << "Deterministic Kernel Weibull: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    if (mortality_tracker != expected_mortality_tracker) {
        cout << "Deterministic Kernel Weibull: mortality tracker (actual, expected):\n"
             << mortality_tracker << "  !=\n"
             << expected_mortality_tracker << "\n";
        return 1;
    }
    return 0;
}

int test_with_log_normal_deterministic_kernel()
{
    Raster<int> infected = {{5, 0, 0}, {0, 5, 0}, {0, 0, 2}};
    Raster<int> mortality_tracker = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    Raster<int> susceptible = {{10, 20, 9}, {14, 15, 0}, {3, 0, 2}};
    Raster<int> total_hosts = susceptible;
    Raster<double> movements = {{0, 0, 1, 1, 2}, {0, 1, 0, 0, 3}, {0, 1, 1, 0, 2}};
    Raster<double> temperature = {{5, 0, 0}, {0, 0, 0}, {0, 0, 2}};
    Raster<double> weather_coefficient = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    std::vector<unsigned> movement_schedule = {1, 1};

    Raster<int> expected_mortality_tracker = {{0, 5, 0}, {5, 0, 0}, {0, 0, 0}};
    Raster<int> expected_infected = {{5, 5, 0}, {5, 5, 0}, {0, 0, 2}};

    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    const std::vector<std::vector<int>> suitable_cell = {
        {0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};

    bool weather = false;
    double reproductive_rate = 2;
    bool generate_stochasticity = false;
    bool establishment_stochasticity = false;
    bool movement_stochasticity = false;
    // We want everything to establish.
    double establishment_probability = 1;
    Simulation<Raster<int>, Raster<double>> s2(
        42,
        infected.rows(),
        infected.cols(),
        model_type_from_string("SI"),
        0,
        generate_stochasticity,
        establishment_stochasticity,
        movement_stochasticity);
    s2.generate(
        dispersers,
        infected,
        weather,
        weather_coefficient,
        reproductive_rate,
        suitable_cell);
    auto expected_dispersers = reproductive_rate * infected;
    if (dispersers != expected_dispersers) {
        cout << "Deterministic Kernel LogNormal: dispersers (actual, expected):\n"
             << dispersers << "  !=\n"
             << expected_dispersers << "\n";
        return 1;
    }
    DeterministicDispersalKernel<Raster<int>> deterministicKernel(
        DispersalKernelType::LogNormal, dispersers, 0.99, 30, 30, 0.1);

    s2.disperse(
        dispersers,
        susceptible,
        infected,
        mortality_tracker,
        total_hosts,
        outside_dispersers,
        weather,
        weather_coefficient,
        deterministicKernel,
        suitable_cell,
        establishment_probability);
    if (outside_dispersers.size() != 8) {
        cout << "Deterministic Kernel LogNormal: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be 8\n";
        return 1;
    }
    if (infected != expected_infected) {
        cout << "Deterministic Kernel LogNormal: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    if (mortality_tracker != expected_mortality_tracker) {
        cout
            << "Deterministic Kernel LogNormal: mortality tracker (actual, expected):\n"
            << mortality_tracker << "  !=\n"
            << expected_mortality_tracker << "\n";
        return 1;
    }
    return 0;
}

int test_with_normal_deterministic_kernel()
{
    Raster<int> infected = {{5, 0, 0}, {0, 5, 0}, {0, 0, 2}};
    Raster<int> mortality_tracker = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    Raster<int> susceptible = {{10, 20, 9}, {14, 15, 0}, {3, 0, 2}};
    Raster<int> total_hosts = susceptible;
    Raster<double> movements = {{0, 0, 1, 1, 2}, {0, 1, 0, 0, 3}, {0, 1, 1, 0, 2}};
    Raster<double> temperature = {{5, 0, 0}, {0, 0, 0}, {0, 0, 2}};
    Raster<double> weather_coefficient = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    std::vector<unsigned> movement_schedule = {1, 1};

    Raster<int> expected_mortality_tracker = {{10, 0, 0}, {0, 10, 0}, {0, 0, 2}};
    Raster<int> expected_infected = {{15, 0, 0}, {0, 15, 0}, {0, 0, 4}};

    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    const std::vector<std::vector<int>> suitable_cell = {
        {0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};

    bool weather = false;
    double reproductive_rate = 2;
    bool generate_stochasticity = false;
    bool establishment_stochasticity = false;
    bool movement_stochasticity = false;
    // We want everything to establish.
    double establishment_probability = 1;
    Simulation<Raster<int>, Raster<double>> s2(
        42,
        infected.rows(),
        infected.cols(),
        model_type_from_string("SI"),
        0,
        generate_stochasticity,
        establishment_stochasticity,
        movement_stochasticity);
    s2.generate(
        dispersers,
        infected,
        weather,
        weather_coefficient,
        reproductive_rate,
        suitable_cell);
    auto expected_dispersers = reproductive_rate * infected;
    if (dispersers != expected_dispersers) {
        cout << "Deterministic Kernel Normal: dispersers (actual, expected):\n"
             << dispersers << "  !=\n"
             << expected_dispersers << "\n";
        return 1;
    }
    DeterministicDispersalKernel<Raster<int>> deterministicKernel(
        DispersalKernelType::Normal, dispersers, 0.99, 30, 30, 1.0);

    s2.disperse(
        dispersers,
        susceptible,
        infected,
        mortality_tracker,
        total_hosts,
        outside_dispersers,
        weather,
        weather_coefficient,
        deterministicKernel,
        suitable_cell,
        establishment_probability);
    if (outside_dispersers.size() != 0) {
        cout << "Deterministic Kernel Normal: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be 0\n";
        return 1;
    }
    if (infected != expected_infected) {
        cout << "Deterministic Kernel Normal: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    if (mortality_tracker != expected_mortality_tracker) {
        cout << "Deterministic Kernel Normal: mortality tracker (actual, expected):\n"
             << mortality_tracker << "  !=\n"
             << expected_mortality_tracker << "\n";
        return 1;
    }
    return 0;
}

int test_with_hyperbolic_secant_deterministic_kernel()
{
    Raster<int> infected = {{5, 0, 0}, {0, 5, 0}, {0, 0, 2}};
    Raster<int> mortality_tracker = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    Raster<int> susceptible = {{10, 20, 9}, {14, 15, 0}, {3, 0, 2}};
    Raster<int> total_hosts = susceptible;
    Raster<double> movements = {{0, 0, 1, 1, 2}, {0, 1, 0, 0, 3}, {0, 1, 1, 0, 2}};
    Raster<double> temperature = {{5, 0, 0}, {0, 0, 0}, {0, 0, 2}};
    Raster<double> weather_coefficient = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    std::vector<unsigned> movement_schedule = {1, 1};

    Raster<int> expected_mortality_tracker = {{10, 0, 0}, {0, 10, 0}, {0, 0, 2}};
    Raster<int> expected_infected = {{15, 0, 0}, {0, 15, 0}, {0, 0, 4}};

    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    const std::vector<std::vector<int>> suitable_cell = {
        {0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};

    bool weather = false;
    double reproductive_rate = 2;
    bool generate_stochasticity = false;
    bool establishment_stochasticity = false;
    bool movement_stochasticity = false;
    // We want everything to establish.
    double establishment_probability = 1;
    Simulation<Raster<int>, Raster<double>> s2(
        42,
        infected.rows(),
        infected.cols(),
        model_type_from_string("SI"),
        0,
        generate_stochasticity,
        establishment_stochasticity,
        movement_stochasticity);
    s2.generate(
        dispersers,
        infected,
        weather,
        weather_coefficient,
        reproductive_rate,
        suitable_cell);
    auto expected_dispersers = reproductive_rate * infected;
    if (dispersers != expected_dispersers) {
        cout
            << "Deterministic Kernel HyperbolicSecant: dispersers (actual, expected):\n"
            << dispersers << "  !=\n"
            << expected_dispersers << "\n";
        return 1;
    }
    DeterministicDispersalKernel<Raster<int>> deterministicKernel(
        DispersalKernelType::HyperbolicSecant, dispersers, 0.99, 30, 30, 1.0);

    s2.disperse(
        dispersers,
        susceptible,
        infected,
        mortality_tracker,
        total_hosts,
        outside_dispersers,
        weather,
        weather_coefficient,
        deterministicKernel,
        suitable_cell,
        establishment_probability);
    if (outside_dispersers.size() != 0) {
        cout << "Deterministic Kernel HyperbolicSecant: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be 0\n";
        return 1;
    }
    if (infected != expected_infected) {
        cout << "Deterministic Kernel HyperbolicSecant: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    if (mortality_tracker != expected_mortality_tracker) {
        cout
            << "Deterministic Kernel HyperbolicSecant: mortality tracker (actual, expected):\n"
            << mortality_tracker << "  !=\n"
            << expected_mortality_tracker << "\n";
        return 1;
    }
    return 0;
}

int test_with_power_law_deterministic_kernel()
{
    Raster<int> infected = {{5, 0, 0}, {0, 5, 0}, {0, 0, 2}};
    Raster<int> mortality_tracker = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    Raster<int> susceptible = {{10, 20, 9}, {14, 15, 0}, {3, 0, 2}};
    Raster<int> total_hosts = susceptible;
    Raster<double> movements = {{0, 0, 1, 1, 2}, {0, 1, 0, 0, 3}, {0, 1, 1, 0, 2}};
    Raster<double> temperature = {{5, 0, 0}, {0, 0, 0}, {0, 0, 2}};
    Raster<double> weather_coefficient = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    std::vector<unsigned> movement_schedule = {1, 1};

    Raster<int> expected_mortality_tracker = {{9, 1, 0}, {0, 9, 0}, {0, 0, 2}};
    Raster<int> expected_infected = {{14, 1, 0}, {0, 14, 0}, {0, 0, 4}};

    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    const std::vector<std::vector<int>> suitable_cell = {
        {0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};

    bool weather = false;
    double reproductive_rate = 2;
    bool generate_stochasticity = false;
    bool establishment_stochasticity = false;
    bool movement_stochasticity = false;
    // We want everything to establish.
    double establishment_probability = 1;
    Simulation<Raster<int>, Raster<double>> s2(
        42,
        infected.rows(),
        infected.cols(),
        model_type_from_string("SI"),
        0,
        generate_stochasticity,
        establishment_stochasticity,
        movement_stochasticity);
    s2.generate(
        dispersers,
        infected,
        weather,
        weather_coefficient,
        reproductive_rate,
        suitable_cell);
    auto expected_dispersers = reproductive_rate * infected;
    if (dispersers != expected_dispersers) {
        cout << "Deterministic Kernel PowerLaw: dispersers (actual, expected):\n"
             << dispersers << "  !=\n"
             << expected_dispersers << "\n";
        return 1;
    }
    DeterministicDispersalKernel<Raster<int>> deterministicKernel(
        DispersalKernelType::PowerLaw, dispersers, 0.99, 30, 30, 1.5, 2);

    s2.disperse(
        dispersers,
        susceptible,
        infected,
        mortality_tracker,
        total_hosts,
        outside_dispersers,
        weather,
        weather_coefficient,
        deterministicKernel,
        suitable_cell,
        establishment_probability);
    if (outside_dispersers.size() != 1) {
        cout << "Deterministic Kernel PowerLaw: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be 0\n";
        return 1;
    }
    if (infected != expected_infected) {
        cout << "Deterministic Kernel PowerLaw: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    if (mortality_tracker != expected_mortality_tracker) {
        cout << "Deterministic Kernel PowerLaw: mortality tracker (actual, expected):\n"
             << mortality_tracker << "  !=\n"
             << expected_mortality_tracker << "\n";
        return 1;
    }
    return 0;
}
int test_with_logistic_deterministic_kernel()
{
    Raster<int> infected = {{5, 0, 0}, {0, 5, 0}, {0, 0, 2}};
    Raster<int> mortality_tracker = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    Raster<int> susceptible = {{10, 20, 9}, {14, 15, 0}, {3, 0, 2}};
    Raster<int> total_hosts = susceptible;
    Raster<double> movements = {{0, 0, 1, 1, 2}, {0, 1, 0, 0, 3}, {0, 1, 1, 0, 2}};
    Raster<double> temperature = {{5, 0, 0}, {0, 0, 0}, {0, 0, 2}};
    Raster<double> weather_coefficient = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    std::vector<unsigned> movement_schedule = {1, 1};

    Raster<int> expected_mortality_tracker = {{10, 0, 0}, {0, 10, 0}, {0, 0, 2}};
    Raster<int> expected_infected = {{15, 0, 0}, {0, 15, 0}, {0, 0, 4}};

    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    const std::vector<std::vector<int>> suitable_cell = {
        {0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};

    bool weather = false;
    double reproductive_rate = 2;
    bool generate_stochasticity = false;
    bool establishment_stochasticity = false;
    bool movement_stochasticity = false;
    // We want everything to establish.
    double establishment_probability = 1;
    Simulation<Raster<int>, Raster<double>> s2(
        42,
        infected.rows(),
        infected.cols(),
        model_type_from_string("SI"),
        0,
        generate_stochasticity,
        establishment_stochasticity,
        movement_stochasticity);
    s2.generate(
        dispersers,
        infected,
        weather,
        weather_coefficient,
        reproductive_rate,
        suitable_cell);
    auto expected_dispersers = reproductive_rate * infected;
    if (dispersers != expected_dispersers) {
        cout << "Deterministic Kernel Logistic: dispersers (actual, expected):\n"
             << dispersers << "  !=\n"
             << expected_dispersers << "\n";
        return 1;
    }
    DeterministicDispersalKernel<Raster<int>> deterministicKernel(
        DispersalKernelType::Logistic, dispersers, 0.99, 30, 30, 1.5);

    s2.disperse(
        dispersers,
        susceptible,
        infected,
        mortality_tracker,
        total_hosts,
        outside_dispersers,
        weather,
        weather_coefficient,
        deterministicKernel,
        suitable_cell,
        establishment_probability);
    if (outside_dispersers.size() != 0) {
        cout << "Deterministic Kernel Logistic: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be 0\n";
        return 1;
    }
    if (infected != expected_infected) {
        cout << "Deterministic Kernel Logistic: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    if (mortality_tracker != expected_mortality_tracker) {
        cout << "Deterministic Kernel Logistic: mortality tracker (actual, expected):\n"
             << mortality_tracker << "  !=\n"
             << expected_mortality_tracker << "\n";
        return 1;
    }
    return 0;
}
int test_with_gamma_deterministic_kernel()
{
    Raster<int> infected = {{5, 0, 0}, {0, 5, 0}, {0, 0, 2}};
    Raster<int> mortality_tracker = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    Raster<int> susceptible = {{10, 20, 9}, {14, 15, 0}, {3, 0, 2}};
    Raster<int> total_hosts = susceptible;
    Raster<double> movements = {{0, 0, 1, 1, 2}, {0, 1, 0, 0, 3}, {0, 1, 1, 0, 2}};
    Raster<double> temperature = {{5, 0, 0}, {0, 0, 0}, {0, 0, 2}};
    Raster<double> weather_coefficient = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    std::vector<unsigned> movement_schedule = {1, 1};

    Raster<int> expected_mortality_tracker = {{0, 5, 0}, {5, 0, 0}, {0, 0, 0}};
    Raster<int> expected_infected = {{5, 5, 0}, {5, 5, 0}, {0, 0, 2}};

    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    const std::vector<std::vector<int>> suitable_cell = {
        {0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};

    bool weather = false;
    double reproductive_rate = 2;
    bool generate_stochasticity = false;
    bool establishment_stochasticity = false;
    bool movement_stochasticity = false;
    // We want everything to establish.
    double establishment_probability = 1;
    Simulation<Raster<int>, Raster<double>> s2(
        42,
        infected.rows(),
        infected.cols(),
        model_type_from_string("SI"),
        0,
        generate_stochasticity,
        establishment_stochasticity,
        movement_stochasticity);
    s2.generate(
        dispersers,
        infected,
        weather,
        weather_coefficient,
        reproductive_rate,
        suitable_cell);
    auto expected_dispersers = reproductive_rate * infected;
    if (dispersers != expected_dispersers) {
        cout << "Deterministic Kernel Gamma: dispersers (actual, expected):\n"
             << dispersers << "  !=\n"
             << expected_dispersers << "\n";
        return 1;
    }
    DeterministicDispersalKernel<Raster<int>> deterministicKernel(
        DispersalKernelType::Gamma, dispersers, 0.9, 30, 30, 1.5, 0.5);

    s2.disperse(
        dispersers,
        susceptible,
        infected,
        mortality_tracker,
        total_hosts,
        outside_dispersers,
        weather,
        weather_coefficient,
        deterministicKernel,
        suitable_cell,
        establishment_probability);
    if (outside_dispersers.size() != 8) {
        cout << "Deterministic Kernel Gamma: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be 8\n";
        return 1;
    }
    if (infected != expected_infected) {
        cout << "Deterministic Kernel Gamma: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    if (mortality_tracker != expected_mortality_tracker) {
        cout << "Deterministic Kernel Gamma: mortality tracker (actual, expected):\n"
             << mortality_tracker << "  !=\n"
             << expected_mortality_tracker << "\n";
        return 1;
    }
    return 0;
}

int test_with_exponential_power_deterministic_kernel()
{
    Raster<int> infected = {{5, 0, 0}, {0, 5, 0}, {0, 0, 2}};
    Raster<int> mortality_tracker = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    Raster<int> susceptible = {{10, 20, 9}, {14, 15, 0}, {3, 0, 2}};
    Raster<int> total_hosts = susceptible;
    Raster<double> movements = {{0, 0, 1, 1, 2}, {0, 1, 0, 0, 3}, {0, 1, 1, 0, 2}};
    Raster<double> temperature = {{5, 0, 0}, {0, 0, 0}, {0, 0, 2}};
    Raster<double> weather_coefficient = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    std::vector<unsigned> movement_schedule = {1, 1};

    Raster<int> expected_mortality_tracker = {{10, 0, 0}, {0, 10, 0}, {0, 0, 2}};
    Raster<int> expected_infected = {{15, 0, 0}, {0, 15, 0}, {0, 0, 4}};

    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    const std::vector<std::vector<int>> suitable_cell = {
        {0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};

    bool weather = false;
    double reproductive_rate = 2;
    bool generate_stochasticity = false;
    bool establishment_stochasticity = false;
    bool movement_stochasticity = false;
    // We want everything to establish.
    double establishment_probability = 1;
    Simulation<Raster<int>, Raster<double>> s2(
        42,
        infected.rows(),
        infected.cols(),
        model_type_from_string("SI"),
        0,
        generate_stochasticity,
        establishment_stochasticity,
        movement_stochasticity);
    s2.generate(
        dispersers,
        infected,
        weather,
        weather_coefficient,
        reproductive_rate,
        suitable_cell);
    auto expected_dispersers = reproductive_rate * infected;
    if (dispersers != expected_dispersers) {
        cout
            << "Deterministic Kernel Exponential Power: dispersers (actual, expected):\n"
            << dispersers << "  !=\n"
            << expected_dispersers << "\n";
        return 1;
    }
    DeterministicDispersalKernel<Raster<int>> deterministicKernel(
        DispersalKernelType::ExponentialPower, dispersers, 0.9, 30, 30, 1.5, 0.5);

    s2.disperse(
        dispersers,
        susceptible,
        infected,
        mortality_tracker,
        total_hosts,
        outside_dispersers,
        weather,
        weather_coefficient,
        deterministicKernel,
        suitable_cell,
        establishment_probability);
    if (outside_dispersers.size() != 0) {
        cout << "Deterministic Kernel Exponential Power: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be 0\n";
        return 1;
    }
    if (infected != expected_infected) {
        cout << "Deterministic Kernel Exponential Power: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    if (mortality_tracker != expected_mortality_tracker) {
        cout
            << "Deterministic Kernel Exponential Power: mortality tracker (actual, expected):\n"
            << mortality_tracker << "  !=\n"
            << expected_mortality_tracker << "\n";
        return 1;
    }
    return 0;
}

int test_cauchy_distribution_functions()
{
    // testing cauchy pdf & icdf
    double scale = 0.001;  // rounding to thousands place
    CauchyKernel cauchy(1.0);
    double probability = (int)(cauchy.pdf(5) / scale) * scale;
    double probability_ref = 0.012;
    if (probability != probability_ref) {
        cout << "Cauchy Distribution: probability was " << probability
             << " but should be " << probability_ref << "\n";
        return 1;
    }
    double x = (int)(cauchy.icdf(0.98) / scale) * scale;
    double x_ref = 15.894;
    if (x != x_ref) {
        cout << "Cauchy Distribution: x was " << x << " but should be " << x_ref
             << "\n";
        return 1;
    }
    CauchyKernel cauchy1(1.5);
    probability_ref = 0.017;
    probability = (int)(cauchy1.pdf(5) / scale) * scale;
    if (probability != probability_ref) {
        cout << "Cauchy Distribution: probability was " << probability
             << " but should be " << probability_ref << "\n";
        return 1;
    }
    x = (int)(cauchy1.icdf(0.98) / scale) * scale;
    x_ref = 23.841;
    if (x != x_ref) {
        cout << "Cauchy Distribution: x was " << x << " but should be " << x_ref
             << "\n";
        return 1;
    }
    return 0;
}

int test_exponential_distribution_functions()
{
    double scale = 0.001;  // rounding to thousands place
    // testing exponential pdf & icdf
    ExponentialKernel exponential(1.0);
    double probability = (int)(exponential.pdf(1) / scale) * scale;
    double probability_ref = 0.367;
    if (probability != probability_ref) {
        cout << "Exponential Distribution: probability was " << probability
             << " but should be " << probability_ref << "\n";
        return 1;
    }
    double x = (int)(exponential.icdf(0.98) / scale) * scale;
    double x_ref = 3.912;
    if (x != x_ref) {
        cout << "Exponential Distribution: x was " << x << " but should be " << x_ref
             << "\n";
        return 1;
    }
    ExponentialKernel exponential2(1.5);
    probability = (int)(exponential2.pdf(1) / scale) * scale;
    probability_ref = 0.342;
    if (probability != probability_ref) {
        cout << "Exponential Distribution: probability was " << probability
             << " but should be " << probability_ref << "\n";
        return 1;
    }
    x = (int)(exponential2.icdf(0.98) / scale) * scale;
    x_ref = 5.868;
    if (x != x_ref) {
        cout << "Exponential Distribution: x was " << x << " but should be " << x_ref
             << "\n";
        return 1;
    }
    return 0;
}

int test_gamma_distribution_functions()
{
    cout << "Gamma\n";
    double x = 0.1;
    for (double a = 0.1; a < 30; a += 5) {
        for (double t = 0.1; t < 30; t += 5) {
            x = 0.1;
            for (double icdf_x = 0.1; icdf_x < 1; icdf_x += 0.1) {
                // testing gamma pdf & icdf
                cout << "alpha = " << a << " theta = " << t << " x = " << x << '\n';
                GammaKernel gamma(a, t);
                double probability = gamma.pdf(x);
                cout << " pdf = " << probability << '\n';
                double icdf = -11;
                try {
                    icdf = gamma.icdf(icdf_x);
                }
                catch (const std::invalid_argument& ia) {
                    std::cerr << "Invalid argument: " << ia.what() << "x = " << x
                              << " alpha = " << a << " theta = " << t << '\n';
                }
                cout << " icdf = " << icdf << "\n";
                x += 5;
            }
        }
    }
    return 1;
}

int test_exponential_power_distribution_functions()
{
    cout << "Exponential Power\n";
    for (double a = 0.1; a < 30; a += 5) {
        for (double b = 0.5; b < 5; b += 1) {
            for (double x = 0.1; x < 1; x += 0.1) {
                // testing gamma pdf & icdf
                ExponentialPowerKernel ep(a, b);
                cout << "alpha = " << a << " beta = " << b << " x = " << x;
                double probability = ep.pdf(x);
                cout << " pdf = " << probability;
                double icdf = ep.icdf(x);
                cout << " icdf = " << icdf << "\n";
            }
        }
    }
    return 1;
}

int test_log_normal_distribution_functions()
{
    cout << "Log Normal\n";
    double x = 0;
    for (double s = 0.1; s < 10; s += 1) {
        x = 0;
        for (double icdf_x = 0.1; icdf_x < 1; icdf_x += 0.1) {
            // testing gamma pdf & icdf
            LogNormalKernel ln(s);
            cout << "sigma = " << s << " x = " << x;
            double probability = ln.pdf(x);
            cout << " pdf = " << probability << " icdf x = " << icdf_x;
            double icdf = ln.icdf(icdf_x);
            cout << " icdf = " << icdf << "\n";
            x += 0.25;
        }
    }
    return 1;
}

int main()
{
    int ret = 0;

    ret += test_with_exponential_deterministic_kernel();
    ret += test_with_cauchy_deterministic_kernel();
    ret += test_cauchy_distribution_functions();
    ret += test_exponential_distribution_functions();
    ret += test_with_weibull_deterministic_kernel();
    ret += test_with_log_normal_deterministic_kernel();
    ret += test_with_normal_deterministic_kernel();
    ret += test_with_hyperbolic_secant_deterministic_kernel();
    ret += test_with_power_law_deterministic_kernel();
    ret += test_with_logistic_deterministic_kernel();
    ret += test_with_gamma_deterministic_kernel();
    ret += test_with_exponential_power_deterministic_kernel();
    // ret += test_gamma_distribution_functions();
    // ret += test_exponential_power_distribution_functions();
    // ret += test_log_normal_distribution_functions();

    std::cout << "Test deterministic number of errors: " << ret << std::endl;
    return ret;
}
#endif  // POPS_TEST
