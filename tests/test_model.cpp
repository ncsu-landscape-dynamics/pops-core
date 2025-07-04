#ifdef POPS_TEST

/*
 * Tests for the PoPS Model class.
 *
 * Copyright (C) 2020 by the authors.
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

#include <vector>

#include <pops/model.hpp>
#include <pops/spread_rate.hpp>

using namespace pops;
using std::cout;

int test_with_reduced_stochasticity()
{
    Raster<int> infected = {{5, 0}, {0, 0}};
    Raster<int> susceptible = {{10, 20}, {14, 15}};
    Raster<int> total_hosts = {{15, 20}, {14, 15}};
    Raster<int> total_populations = {{20, 20}, {20, 20}};
    Raster<int> zeros(infected.rows(), infected.cols(), 0);

    Raster<int> expected_mortality_tracker = {{0, 10}, {0, 0}};
    auto expected_infected = expected_mortality_tracker + infected;

    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};

    Raster<int> dispersers(infected.rows(), infected.cols());
    Raster<int> established_dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;
    Config config;
    config.weather = false;
    config.reproductive_rate = 2;
    config.generate_stochasticity = false;
    config.establishment_stochasticity = false;
    // We want everything to establish.
    config.establishment_probability = 1;
    config.natural_kernel_type = "deterministic neighbor";
    config.natural_direction = "E";
    config.use_anthropogenic_kernel = false;
    config.random_seed = 42;
    config.rows = infected.rows();
    config.cols = infected.cols();
    config.model_type = "SI";
    config.latency_period_steps = 0;
    config.use_lethal_temperature = false;
    config.use_survival_rate = false;
    config.use_quarantine = true;
    config.quarantine_frequency = "year";
    config.quarantine_frequency_n = 1;
    config.use_spreadrates = false;
    config.natural_scale = 0.9;
    config.anthro_scale = 0.9;
    config.set_date_start(2020, 1, 1);
    config.set_date_end(2021, 12, 31);
    config.set_step_unit(StepUnit::Month);
    config.set_step_num_units(1);
    config.use_mortality = true;
    config.mortality_frequency = "year";
    config.mortality_frequency_n = 1;
    config.use_treatments = false;
    config.ew_res = 1;
    config.ns_res = 1;
    config.create_schedules();

    unsigned num_mortality_steps = 1;
    std::vector<Raster<int>> mortality_tracker(
        num_mortality_steps, Raster<int>(infected.rows(), infected.cols(), 0));

    //    int exposed_size = 0;
    //    if (config.latency_period_steps)
    //        exposed_size = config.latency_period_steps + 1;
    //    std::vector<Raster<int>> exposed(
    //                exposed_size,
    //                Raster<int>(infected.rows(), infected.cols(), 0));
    Raster<int> died(infected.rows(), infected.cols(), 0);
    Raster<int> total_exposed(infected.rows(), infected.cols(), 0);
    std::vector<Raster<int>> empty_integer;
    std::vector<Raster<double>> empty_float;

    unsigned quarantine_num_steps =
        get_number_of_scheduled_actions(config.quarantine_schedule());
    QuarantineEscapeAction<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, quarantine_num_steps);

    auto expected_dispersers = config.reproductive_rate * infected;
    auto expected_established_dispersers = config.reproductive_rate * infected;
    std::vector<std::vector<int>> movements = {
        {0, 0, 1, 1, 2}, {0, 1, 0, 0, 3}, {0, 1, 1, 0, 2}};

    int step = 0;

    Model<Raster<int>, Raster<double>, Raster<double>::IndexType> model(config);
    model.run_step(
        step++,
        infected,
        susceptible,
        total_populations,
        total_hosts,
        dispersers,
        established_dispersers,
        total_exposed,
        empty_integer,
        mortality_tracker,
        died,
        empty_float,
        empty_float,
        zeros,
        outside_dispersers,
        quarantine,
        zeros,
        movements,
        Network<int>::null_network(),
        suitable_cells);

    int ret = 0;
    if (dispersers != expected_dispersers) {
        cout << "reduced_stochasticity: dispersers (actual, expected):\n"
             << dispersers << "  !=\n"
             << expected_dispersers << "\n";
        ++ret;
    }
    if (established_dispersers != expected_established_dispersers) {
        cout << "reduced_stochasticity: established dispersers (actual, expected):\n"
             << established_dispersers << "  !=\n"
             << expected_established_dispersers << "\n";
        ++ret;
    }
    if (!outside_dispersers.empty()) {
        cout << "reduced_stochasticity: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be none\n";
        ++ret;
    }
    if (infected != expected_infected) {
        cout << "reduced_stochasticity: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        ++ret;
    }
    if (mortality_tracker[0] != expected_mortality_tracker) {
        cout << "reduced_stochasticity: mortality tracker (actual, expected):\n"
             << mortality_tracker[0] << "  !=\n"
             << expected_mortality_tracker << "\n";
        ++ret;
    }
    return ret;
}
int test_deterministic()
{
    Raster<int> infected = {{5, 0, 0}, {0, 5, 0}, {0, 0, 2}};
    Raster<int> susceptible = {{10, 20, 9}, {14, 15, 0}, {3, 0, 2}};
    Raster<int> total_hosts = susceptible + infected;
    Raster<int> total_populations = {{20, 20, 20}, {20, 20, 20}, {20, 20, 20}};
    Raster<int> zeros(infected.rows(), infected.cols(), 0);

    Raster<int> expected_mortality_tracker = {{10, 0, 0}, {0, 10, 0}, {0, 0, 2}};
    Raster<int> expected_infected = {{15, 0, 0}, {0, 15, 0}, {0, 0, 4}};

    Raster<int> dispersers(infected.rows(), infected.cols());
    Raster<int> established_dispersers(infected.rows(), infected.cols());

    std::vector<std::tuple<int, int>> outside_dispersers;

    std::vector<std::vector<int>> suitable_cells = {
        {0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};

    Config config;
    config.weather = false;
    config.reproductive_rate = 2;
    config.generate_stochasticity = false;
    config.establishment_stochasticity = false;
    config.movement_stochasticity = false;
    std::vector<std::vector<int>> movements = {
        {0, 0, 1, 1, 2}, {0, 1, 0, 0, 3}, {0, 1, 1, 0, 2}};
    config.movement_schedule = {1, 1};

    // We want everything to establish.
    config.establishment_probability = 1;
    config.natural_kernel_type = "cauchy";
    config.natural_direction = "none";
    config.natural_scale = 0.9;
    config.anthro_scale = 0.9;
    config.dispersal_percentage = 0.9;
    config.natural_kappa = 0;
    config.anthro_kappa = 0;

    config.use_anthropogenic_kernel = false;
    config.random_seed = 42;
    config.rows = infected.rows();
    config.cols = infected.cols();
    config.model_type = "SI";
    config.latency_period_steps = 0;
    config.use_lethal_temperature = false;
    config.use_survival_rate = false;
    config.use_quarantine = false;
    config.use_spreadrates = false;

    config.set_date_start(2020, 1, 1);
    config.set_date_end(2021, 12, 31);
    config.set_step_unit(StepUnit::Month);
    config.set_step_num_units(1);
    config.use_mortality = true;
    config.mortality_frequency = "year";
    config.mortality_frequency_n = 1;
    config.use_treatments = false;
    config.ew_res = 30;
    config.ns_res = 30;
    config.create_schedules();

    config.dispersal_stochasticity = false;

    unsigned num_mortality_steps = 1;
    std::vector<Raster<int>> mortality_tracker(
        num_mortality_steps, Raster<int>(infected.rows(), infected.cols(), 0));

    Raster<int> died(infected.rows(), infected.cols(), 0);
    Raster<int> total_exposed(infected.rows(), infected.cols(), 0);
    std::vector<Raster<int>> empty_integer;
    std::vector<Raster<double>> empty_floats;
    QuarantineEscapeAction<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, 0);

    auto expected_dispersers = config.reproductive_rate * infected;
    auto expected_established_dispersers = config.reproductive_rate * infected;

    // Limit established dispersers by number of available hosts.
    for (int row = 0; row < susceptible.rows(); ++row) {
        for (int col = 0; col < susceptible.cols(); ++col) {
            if (expected_established_dispersers(row, col) > susceptible(row, col)) {
                expected_established_dispersers(row, col) = susceptible(row, col);
            }
        }
    }

    int step = 0;

    Model<Raster<int>, Raster<double>, Raster<double>::IndexType> model(config);
    model.run_step(
        step++,
        infected,
        susceptible,
        total_populations,
        total_hosts,
        dispersers,
        established_dispersers,
        total_exposed,
        empty_integer,
        mortality_tracker,
        died,
        empty_floats,
        empty_floats,
        zeros,
        outside_dispersers,
        quarantine,
        zeros,
        movements,
        Network<int>::null_network(),
        suitable_cells);

    if (dispersers != expected_dispersers) {
        cout << "deterministic: dispersers (actual, expected):\n"
             << dispersers << "  !=\n"
             << expected_dispersers << "\n";
        return 1;
    }
    if (established_dispersers != expected_established_dispersers) {
        cout << "deterministic: established dispersers (actual, expected):\n"
             << established_dispersers << "  !=\n"
             << expected_established_dispersers << "\n";
        return 1;
    }
    if (!outside_dispersers.empty()) {
        cout << "deterministic: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be none\n";
        return 1;
    }
    if (infected != expected_infected) {
        cout << "deterministic: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    if (mortality_tracker[0] != expected_mortality_tracker) {
        cout << "deterministic: mortality tracker (actual, expected):\n"
             << mortality_tracker[0] << "  !=\n"
             << expected_mortality_tracker << "\n";
        return 1;
    }
    return 0;
}
int test_deterministic_exponential()
{
    Raster<int> infected = {{5, 0, 0}, {0, 5, 0}, {0, 0, 2}};
    Raster<int> susceptible = {{10, 20, 9}, {14, 15, 0}, {3, 0, 2}};
    Raster<int> total_hosts = susceptible + infected;
    Raster<int> total_populations = {{20, 20, 20}, {20, 20, 20}, {20, 20, 20}};
    Raster<int> zeros(infected.rows(), infected.cols(), 0);

    Raster<int> expected_mortality_tracker = {{10, 0, 0}, {0, 10, 0}, {0, 0, 2}};
    Raster<int> expected_infected = {{15, 0, 0}, {0, 15, 0}, {0, 0, 4}};

    Raster<int> dispersers(infected.rows(), infected.cols());
    Raster<int> established_dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    std::vector<std::vector<int>> suitable_cells = {
        {0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};

    Config config;
    config.weather = false;
    config.reproductive_rate = 2;
    config.generate_stochasticity = false;
    config.establishment_stochasticity = false;
    config.movement_stochasticity = false;
    std::vector<std::vector<int>> movements = {
        {0, 0, 1, 1, 2}, {0, 1, 0, 0, 3}, {0, 1, 1, 0, 2}};
    config.movement_schedule = {1, 1};

    // We want everything to establish.
    config.establishment_probability = 1;
    config.natural_kernel_type = "exponential";
    config.natural_direction = "none";
    config.natural_scale = 1;
    config.anthro_scale = 1;
    config.natural_kappa = 0;
    config.anthro_kappa = 0;
    config.dispersal_percentage = 0.99;
    config.use_anthropogenic_kernel = false;
    config.random_seed = 42;
    config.rows = infected.rows();
    config.cols = infected.cols();
    config.model_type = "SI";
    config.latency_period_steps = 0;
    config.use_lethal_temperature = false;
    config.use_survival_rate = false;
    config.use_quarantine = false;
    config.use_spreadrates = false;

    config.set_date_start(2020, 1, 1);
    config.set_date_end(2021, 12, 31);
    config.set_step_unit(StepUnit::Month);
    config.set_step_num_units(1);
    config.use_mortality = true;
    config.mortality_frequency = "year";
    config.mortality_frequency_n = 1;
    config.use_treatments = false;
    config.ew_res = 30;
    config.ns_res = 30;
    config.create_schedules();

    config.dispersal_stochasticity = false;

    unsigned num_mortality_steps = 1;
    std::vector<Raster<int>> mortality_tracker(
        num_mortality_steps, Raster<int>(infected.rows(), infected.cols(), 0));

    Raster<int> died(infected.rows(), infected.cols(), 0);
    Raster<int> total_exposed(infected.rows(), infected.cols(), 0);
    std::vector<Raster<int>> empty_integer;
    std::vector<Raster<double>> empty_floats;
    QuarantineEscapeAction<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, 0);

    auto expected_dispersers = config.reproductive_rate * infected;
    auto expected_established_dispersers = config.reproductive_rate * infected;

    // Limit established dispersers by number of available hosts.
    for (int row = 0; row < susceptible.rows(); ++row) {
        for (int col = 0; col < susceptible.cols(); ++col) {
            if (expected_established_dispersers(row, col) > susceptible(row, col)) {
                expected_established_dispersers(row, col) = susceptible(row, col);
            }
        }
    }

    int step = 0;

    Model<Raster<int>, Raster<double>, Raster<double>::IndexType> model(config);
    model.run_step(
        step++,
        infected,
        susceptible,
        total_populations,
        total_hosts,
        dispersers,
        established_dispersers,
        total_exposed,
        empty_integer,
        mortality_tracker,
        died,
        empty_floats,
        empty_floats,
        zeros,
        outside_dispersers,
        quarantine,
        zeros,
        movements,
        Network<int>::null_network(),
        suitable_cells);

    if (dispersers != expected_dispersers) {
        cout << "deterministic exponential: dispersers (actual, expected):\n"
             << dispersers << "  !=\n"
             << expected_dispersers << "\n";
        return 1;
    }
    if (established_dispersers != expected_established_dispersers) {
        cout
            << "deterministic exponential: established dispersers (actual, expected):\n"
            << established_dispersers << "  !=\n"
            << expected_established_dispersers << "\n";
        return 1;
    }
    if (!outside_dispersers.empty()) {
        cout << "deterministic exponential: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be none\n";
        return 1;
    }
    if (infected != expected_infected) {
        cout << "deterministic exponential: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    if (mortality_tracker[0] != expected_mortality_tracker) {
        cout << "deterministic exponential: mortality tracker (actual, expected):\n"
             << mortality_tracker[0] << "  !=\n"
             << expected_mortality_tracker << "\n";
        return 1;
    }
    return 0;
}

int test_model_sei_deterministic()
{
    Raster<int> infected = {{5, 0, 0}, {0, 5, 0}, {0, 0, 2}};
    Raster<int> susceptible = {{95, 100, 100}, {100, 95, 100}, {100, 0, 98}};
    Raster<int> total_hosts = susceptible + infected;
    Raster<int> total_populations = {{100, 100, 100}, {100, 100, 100}, {100, 100, 100}};
    Raster<int> zeros(infected.rows(), infected.cols(), 0);

    Raster<int> dispersers(infected.rows(), infected.cols());
    Raster<int> established_dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    std::vector<std::vector<int>> suitable_cells = {
        {0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};

    Config config;
    config.reproductive_rate = 1;
    config.generate_stochasticity = false;
    config.establishment_stochasticity = false;
    config.movement_stochasticity = false;

    // We want everything to establish.
    config.establishment_probability = 1;
    config.natural_kernel_type = "cauchy";
    config.natural_direction = "none";
    config.natural_scale = 0.9;
    config.anthro_scale = 0.9;
    config.dispersal_percentage = 0.9;
    config.natural_kappa = 0;
    config.anthro_kappa = 0;

    config.use_anthropogenic_kernel = false;
    config.random_seed = 42;
    config.rows = infected.rows();
    config.cols = infected.cols();
    config.model_type = "SEI";
    config.latency_period_steps = 11;
    config.use_lethal_temperature = false;
    config.use_survival_rate = false;
    config.use_quarantine = false;
    config.use_spreadrates = false;

    config.set_date_start(2020, 1, 1);
    config.set_date_end(2020, 12, 31);
    config.set_step_unit(StepUnit::Month);
    config.set_step_num_units(1);
    config.use_mortality = false;
    config.mortality_frequency = "year";
    config.mortality_frequency_n = 1;
    config.use_treatments = false;
    config.ew_res = 30;
    config.ns_res = 30;
    config.create_schedules();

    config.dispersal_stochasticity = false;

    std::vector<std::vector<int>> movements;

    unsigned num_mortality_steps = 1;
    std::vector<Raster<int>> mortality_tracker(
        num_mortality_steps, Raster<int>(infected.rows(), infected.cols(), 0));

    Raster<int> died(infected.rows(), infected.cols(), 0);
    Raster<int> total_exposed(infected.rows(), infected.cols(), 0);
    int exposed_size = 0;
    if (config.latency_period_steps)
        exposed_size = config.latency_period_steps + 1;
    std::vector<Raster<int>> exposed(
        exposed_size, Raster<int>(infected.rows(), infected.cols(), 0));
    std::vector<Raster<double>> empty_float;
    QuarantineEscapeAction<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, 0);

    // There should be still the original number of infected when dispersers are
    // created.
    auto expected_dispersers = config.reproductive_rate * infected;
    auto expected_established_dispersers = config.reproductive_rate * infected;
    // One E to I transition should happen.
    auto expected_infected = config.reproductive_rate * infected + infected;

    Model<Raster<int>, Raster<double>, Raster<double>::IndexType> model(config);
    for (unsigned int step = 0; step < config.scheduler().get_num_steps(); ++step) {
        model.run_step(
            step,
            infected,
            susceptible,
            total_populations,
            total_hosts,
            dispersers,
            established_dispersers,
            total_exposed,
            exposed,
            mortality_tracker,
            died,
            empty_float,
            empty_float,
            zeros,
            outside_dispersers,
            quarantine,
            zeros,
            movements,
            Network<int>::null_network(),
            suitable_cells);
    }
    if (dispersers != expected_dispersers) {
        cout << "sei_deterministic: dispersers (actual, expected):\n"
             << dispersers << "  !=\n"
             << expected_dispersers << "\n";
        return 1;
    }
    if (established_dispersers != expected_established_dispersers) {
        cout
            << "sei_deterministic exponential: established dispersers (actual, expected):\n"
            << established_dispersers << "  !=\n"
            << expected_established_dispersers << "\n";
        return 1;
    }
    if (!outside_dispersers.empty()) {
        cout << "sei_deterministic: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be none\n";
        return 1;
    }
    if (infected != expected_infected) {
        cout << "sei_deterministic: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    return 0;
}

int test_model_sei_deterministic_with_treatments()
{
    int ret = 0;
    Raster<int> infected = {{5, 0, 0}, {0, 10, 0}, {0, 0, 2}};
    Raster<int> susceptible = {{95, 100, 100}, {100, 95, 100}, {100, 0, 98}};
    Raster<int> total_hosts = susceptible + infected;
    Raster<int> total_populations = {{100, 100, 100}, {100, 100, 100}, {100, 100, 100}};
    Raster<int> zeros(infected.rows(), infected.cols(), 0);

    Raster<int> dispersers(infected.rows(), infected.cols());
    Raster<int> established_dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    std::vector<std::vector<int>> suitable_cells = {
        {0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};

    Config config;
    config.reproductive_rate = 1;
    config.generate_stochasticity = false;
    config.establishment_stochasticity = false;
    config.movement_stochasticity = false;

    // We want everything to establish.
    config.establishment_probability = 1;
    config.natural_kernel_type = "cauchy";
    config.natural_direction = "none";
    config.natural_scale = 0.9;
    config.anthro_scale = 0.9;
    config.dispersal_percentage = 0.9;
    config.natural_kappa = 0;
    config.anthro_kappa = 0;

    config.use_anthropogenic_kernel = false;
    config.random_seed = 42;
    config.rows = infected.rows();
    config.cols = infected.cols();
    config.model_type = "SEI";
    config.latency_period_steps = 11;
    config.use_lethal_temperature = false;
    config.use_survival_rate = false;
    config.use_quarantine = false;
    config.use_spreadrates = false;

    config.set_date_start(2020, 1, 1);
    config.set_date_end(2020, 12, 31);
    config.set_step_unit(StepUnit::Month);
    config.set_step_num_units(1);
    config.use_mortality = true;
    config.mortality_frequency = "year";
    config.mortality_frequency_n = 1;
    config.use_treatments = true;
    config.ew_res = 30;
    config.ns_res = 30;
    config.create_schedules();

    config.dispersal_stochasticity = false;

    using TestModel = Model<Raster<int>, Raster<double>, Raster<double>::IndexType>;
    TestModel model{config};

    unsigned num_mortality_steps = 1;
    std::vector<Raster<int>> mortality_tracker(
        num_mortality_steps, Raster<int>(infected.rows(), infected.cols(), 0));

    Raster<int> died(infected.rows(), infected.cols(), 0);
    Raster<int> total_exposed(infected.rows(), infected.cols(), 0);
    Raster<int> resistant(infected.rows(), infected.cols(), 0);
    int exposed_size = 0;
    if (config.latency_period_steps)
        exposed_size = config.latency_period_steps + 1;
    std::vector<Raster<int>> exposed(
        exposed_size, Raster<int>(infected.rows(), infected.cols(), 0));
    std::vector<Raster<double>> empty_floats;
    std::vector<std::vector<int>> movements;

    TestModel::StandardSingleHostPool host_pool(
        config,
        susceptible,
        exposed,
        infected,
        total_exposed,
        resistant,
        mortality_tracker,
        died,
        total_hosts,
        model.environment(),
        suitable_cells);
    std::vector<TestModel::StandardSingleHostPool*> host_pools = {&host_pool};
    TestModel::StandardMultiHostPool multi_host_pool(host_pools, config);
    PestHostTable<TestModel::StandardSingleHostPool> pest_host_table(
        model.environment());
    pest_host_table.add_host_info(
        config.establishment_probability,  // using as host susceptibility
        config.mortality_rate,
        config.mortality_time_lag);
    multi_host_pool.set_pest_host_table(pest_host_table);
    TestModel::StandardPestPool pest_pool{
        dispersers, established_dispersers, outside_dispersers};
    SpreadRateAction<TestModel::StandardMultiHostPool, int> spread_rate(
        multi_host_pool, config.rows, config.cols, config.ew_res, config.ns_res, 0);

    Treatments<TestModel::StandardSingleHostPool, Raster<double>> treatments(
        config.scheduler());
    Raster<double> simple_treatment = {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    treatments.add_treatment(
        simple_treatment, Date(2020, 1, 1), 0, TreatmentApplication::AllInfectedInCell);
    Raster<double> pesticide_treatment = {{0, 0, 0}, {0, 0.5, 0}, {0, 0, 0}};
    treatments.add_treatment(
        pesticide_treatment, Date(2020, 1, 1), 365, TreatmentApplication::Ratio);

    QuarantineEscapeAction<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, 0);

    // One E to I transition should happen.
    auto expected_infected = config.reproductive_rate * infected + infected;
    // Apply treatment to expected results (assuming rate == 1)
    // (modifying int with double is not allowed in Raster, so we have to be explicit)
    // Remove infected
    for (int row = 0; row < expected_infected.rows(); ++row)
        for (int col = 0; col < expected_infected.rows(); ++col)
            expected_infected(row, col) *=
                !static_cast<bool>(simple_treatment(row, col));
    // Reduced number of infected
    // (assuming 1 infection (E to I) step completed, i.e. 1 initial state + 1 step)
    for (int row = 0; row < expected_infected.rows(); ++row)
        for (int col = 0; col < expected_infected.rows(); ++col)
            if (pesticide_treatment(row, col) > 0)
                expected_infected(row, col) = std::lround(
                    2 * pesticide_treatment(row, col) * expected_infected(row, col));
    expected_infected(0, 0) += 5;  // based on what is considered a correct result
    expected_infected(1, 1) -= 5;  // based on what is considered a correct result
    // Values are based on the result which is considered correct.
    Raster<int> expected_dispersers = {{5, 0, 0}, {0, 10, 0}, {0, 0, 2}};

    for (unsigned int step = 0; step < config.scheduler().get_num_steps(); ++step) {
        model.run_step(
            step,
            multi_host_pool,
            pest_pool,
            total_populations,
            treatments,
            empty_floats,
            empty_floats,
            spread_rate,
            quarantine,
            zeros,
            movements,
            Network<int>::null_network());
    }
    if (!outside_dispersers.empty()) {
        cout << "sei_deterministic_with_treatments: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be none\n";
        ++ret;
    }
    if (infected != expected_infected) {
        cout << "sei_deterministic_with_treatments: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        ++ret;
    }
    if (dispersers != expected_dispersers) {
        cout << "sei_deterministic_with_treatments: dispersers (actual, expected):\n"
             << dispersers << "  !=\n"
             << expected_dispersers << "\n";
        ++ret;
    }
    return ret;
}

int main()
{
    int ret = 0;

    ret += test_with_reduced_stochasticity();
    ret += test_deterministic();
    ret += test_deterministic_exponential();
    ret += test_model_sei_deterministic();
    ret += test_model_sei_deterministic_with_treatments();
    std::cout << "Test model number of errors: " << ret << std::endl;

    return ret;
}

#endif  // POPS_TEST
