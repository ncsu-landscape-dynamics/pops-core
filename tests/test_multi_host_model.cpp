#ifdef POPS_TEST

/*
 * Tests for the PoPS Model class multi-host functionality.
 *
 * Copyright (C) 2023 by the authors.
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
#include <pops/competency_table.hpp>
#include <pops/pest_host_table.hpp>

using namespace pops;

/** Test multihost with minimal parameters and one host pool only.
 *
 * Values determined by running the test.
 */
int test_minimal_parameters_one_host()
{
    int ret = 0;

    Config config;
    config.model_type = "SI";
    config.reproductive_rate = 2;
    config.establishment_probability = 1;
    config.random_seed = 42;
    config.natural_scale = 20;
    config.natural_kernel_type = "deterministic-neighbor";
    config.natural_direction = "E";
    config.use_anthropogenic_kernel = false;
    config.anthro_scale = 0.9;
    config.anthro_kappa = 0;
    config.use_spreadrates = false;
    config.use_quarantine = true;
    config.create_schedules();

    using TestModel = Model<Raster<int>, Raster<double>, Raster<double>::IndexType>;
    TestModel model{config};

    int step = 0;

    Raster<int> infected = {{5, 0, 0}, {0, 5, 0}, {0, 0, 2}};
    Raster<int> susceptible = {{10, 20, 9}, {14, 15, 0}, {3, 0, 2}};
    Raster<int> total_populations = {{20, 20, 20}, {20, 20, 20}, {20, 20, 20}};
    Raster<int> total_hosts = susceptible + infected;
    std::vector<std::vector<int>> suitable_cells =
        find_suitable_cells<Raster<int>::IndexType, Raster<int>>(total_hosts);

    Raster<int> dispersers(infected.rows(), infected.cols(), 0);
    Raster<int> established_dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    Raster<int> total_exposed(infected.rows(), infected.cols(), 0);
    Raster<int> died(infected.rows(), infected.cols(), 0);

    Raster<int> empty_integer;
    std::vector<Raster<int>> empty_integers;
    std::vector<Raster<double>> empty_floats;

    unsigned num_mortality_steps = 1;
    std::vector<Raster<int>> mortality_tracker(
        num_mortality_steps, Raster<int>(infected.rows(), infected.cols(), 0));

    std::vector<std::vector<int>> movements = {};

    config.ew_res = 30;
    config.ns_res = 30;
    config.rows = infected.rows();
    config.cols = infected.cols();

    Raster<double> weather = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    TestModel::StandardSingleHostPool host_pool(
        model_type_from_string(config.model_type),
        susceptible,
        empty_integers,
        config.latency_period_steps,
        infected,
        total_exposed,
        empty_integer,
        mortality_tracker,
        died,
        total_hosts,
        model.environment(),
        config.generate_stochasticity,
        config.reproductive_rate,
        config.establishment_stochasticity,
        config.establishment_probability,
        config.rows,
        config.cols,
        suitable_cells);
    std::vector<TestModel::StandardSingleHostPool*> host_pools = {&host_pool};
    TestModel::StandardMultiHostPool multi_host_pool(host_pools, config);
    TestModel::StandardPestPool pest_pool{
        dispersers, established_dispersers, outside_dispersers};
    Treatments<TestModel::StandardSingleHostPool, Raster<double>> treatments(
        config.scheduler());
    SpreadRateAction<TestModel::StandardMultiHostPool, int> spread_rate(
        multi_host_pool, config.rows, config.cols, config.ew_res, config.ns_res, 0);
    Raster<int> zeros(infected.rows(), infected.cols(), 0);
    QuarantineEscapeAction<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, 0);
    model.environment().update_weather_coefficient(weather);
    Raster<int> expected_infected = infected;
    model.run_step(
        step++,
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
    expected_infected += Raster<int>({{0, 7, 0}, {0, 0, 0}, {0, 0, 0}});
    if (infected != expected_infected) {
        std::cerr << "test_minimal_parameters_one_host (step " << step
                  << ") infected (actual, expected):\n"
                  << infected << "  !=\n"
                  << expected_infected << "\n";
        ++ret;
    }
    Raster<int> expected_dispersers = {{8, 0, 0}, {0, 8, 0}, {0, 0, 7}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_minimal_parameters_one_host (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    size_t expected_outside_dispersers_size = 7;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_minimal_parameters_one_host (step " << step
                  << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    model.run_step(
        step++,
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
    expected_infected += Raster<int>({{0, 7, 2}, {0, 0, 0}, {0, 0, 0}});
    if (infected != expected_infected) {
        std::cerr << "test_minimal_parameters_one_host (step " << step
                  << ") infected (actual, expected):\n"
                  << infected << "  !=\n"
                  << expected_infected << "\n";
        ++ret;
    }
    expected_dispersers = {{13, 12, 0}, {0, 7, 0}, {0, 0, 1}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_minimal_parameters_one_host (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    expected_outside_dispersers_size = 8;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_minimal_parameters_one_host (step " << step
                  << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    return ret;
}

/** Test two hosts without competency table.
 *
 * Values determined by running the test.
 */
int test_minimal_parameters_two_hosts()
{
    int ret = 0;

    Config config;
    config.model_type = "SI";
    config.reproductive_rate = 2;
    config.establishment_probability = 1;
    config.random_seed = 42;
    config.natural_scale = 20;
    config.natural_kernel_type = "deterministic-neighbor";
    config.natural_direction = "E";
    config.use_anthropogenic_kernel = false;
    config.anthro_scale = 0.9;
    config.anthro_kappa = 0;
    config.use_spreadrates = false;
    config.use_quarantine = true;
    config.create_schedules();

    using TestModel = Model<Raster<int>, Raster<double>, Raster<double>::IndexType>;
    TestModel model{config};

    int step = 0;

    Raster<int> total_hosts_1 = {{10, 20, 9}, {0, 0, 0}, {3, 50, 2}};
    Raster<int> total_hosts_2 = {{10, 20, 9}, {14, 15, 0}, {0, 0, 0}};
    Raster<int> infected_1 = {{5, 0, 0}, {0, 0, 0}, {0, 10, 2}};
    Raster<int> infected_2 = {{5, 0, 0}, {0, 5, 0}, {0, 0, 0}};
    Raster<int> total_populations = {{2, 0, 1}, {0, 2, 3}, {2, 10, 5}};
    total_populations += total_hosts_1;
    total_populations += total_hosts_2;
    Raster<int> susceptible_1 = total_hosts_1 - infected_1;
    Raster<int> susceptible_2 = total_hosts_2 - infected_2;
    std::vector<std::vector<int>> suitable_cells =
        find_suitable_cells<Raster<int>::IndexType, Raster<int>>(
            {&total_hosts_1, &total_hosts_2});

    Raster<int> dispersers(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> established_dispersers(infected_1.rows(), infected_1.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    Raster<int> total_exposed_1(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> total_exposed_2(infected_2.rows(), infected_2.cols(), 0);
    Raster<int> died_1(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> died_2(infected_2.rows(), infected_2.cols(), 0);

    Raster<int> empty_integer;
    std::vector<Raster<int>> empty_integers;
    std::vector<Raster<double>> empty_floats;

    unsigned num_mortality_steps = 1;
    std::vector<Raster<int>> mortality_tracker_1(
        num_mortality_steps, Raster<int>(infected_1.rows(), infected_1.cols(), 0));
    std::vector<Raster<int>> mortality_tracker_2(
        num_mortality_steps, Raster<int>(infected_2.rows(), infected_2.cols(), 0));

    std::vector<std::vector<int>> movements = {};

    config.ew_res = 30;
    config.ns_res = 30;
    config.rows = infected_1.rows();
    config.cols = infected_1.cols();

    Raster<double> weather = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    TestModel::StandardSingleHostPool host_pool_1(
        model_type_from_string(config.model_type),
        susceptible_1,
        empty_integers,
        config.latency_period_steps,
        infected_1,
        total_exposed_1,
        empty_integer,
        mortality_tracker_1,
        died_1,
        total_hosts_1,
        model.environment(),
        config.generate_stochasticity,
        config.reproductive_rate,
        config.establishment_stochasticity,
        config.establishment_probability,
        config.rows,
        config.cols,
        suitable_cells);
    TestModel::StandardSingleHostPool host_pool_2(
        model_type_from_string(config.model_type),
        susceptible_2,
        empty_integers,
        config.latency_period_steps,
        infected_2,
        total_exposed_2,
        empty_integer,
        mortality_tracker_2,
        died_2,
        total_hosts_2,
        model.environment(),
        config.generate_stochasticity,
        config.reproductive_rate,
        config.establishment_stochasticity,
        config.establishment_probability,
        config.rows,
        config.cols,
        suitable_cells);
    std::vector<TestModel::StandardSingleHostPool*> host_pools = {
        &host_pool_1, &host_pool_2};
    TestModel::StandardMultiHostPool multi_host_pool(host_pools, config);
    TestModel::StandardPestPool pest_pool{
        dispersers, established_dispersers, outside_dispersers};
    Treatments<TestModel::StandardSingleHostPool, Raster<double>> treatments(
        config.scheduler());
    SpreadRateAction<TestModel::StandardMultiHostPool, int> spread_rate(
        multi_host_pool, config.rows, config.cols, config.ew_res, config.ns_res, 0);
    Raster<int> zeros(infected_1.rows(), infected_1.cols(), 0);
    QuarantineEscapeAction<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, 0);

    model.environment().update_weather_coefficient(weather);

    Raster<int> expected_infected_1 = infected_1;
    Raster<int> expected_infected_2 = infected_2;
    model.run_step(
        step++,
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
    expected_infected_1 += Raster<int>({{0, 4, 0}, {0, 0, 0}, {0, 0, 0}});
    if (infected_1 != expected_infected_1) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
                  << ") infected 1 (actual, expected):\n"
                  << infected_1 << "  !=\n"
                  << expected_infected_1 << "\n";
        ++ret;
    }
    expected_infected_2 += Raster<int>({{0, 5, 0}, {0, 0, 0}, {0, 0, 0}});
    if (infected_2 != expected_infected_2) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
                  << ") infected 2 (actual, expected):\n"
                  << infected_2 << "  !=\n"
                  << expected_infected_2 << "\n";
        ++ret;
    }
    Raster<int> expected_dispersers = {{16, 0, 0}, {0, 10, 0}, {0, 22, 5}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    size_t expected_outside_dispersers_size = 5;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
                  << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    model.run_step(
        step++,
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
    expected_infected_1 += Raster<int>({{0, 5, 3}, {0, 0, 0}, {0, 0, 0}});
    if (infected_1 != expected_infected_1) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
                  << ") infected 1 (actual, expected):\n"
                  << infected_1 << "  !=\n"
                  << expected_infected_1 << "\n";
        ++ret;
    }
    expected_infected_2 += Raster<int>({{0, 7, 3}, {0, 0, 0}, {0, 0, 0}});
    if (infected_2 != expected_infected_2) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
                  << ") infected 2 (actual, expected):\n"
                  << infected_2 << "  !=\n"
                  << expected_infected_2 << "\n";
        ++ret;
    }
    expected_dispersers = {{31, 20, 0}, {0, 13, 0}, {0, 13, 6}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    expected_outside_dispersers_size = 11;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
                  << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    if (ret) {
        std::cerr << "Next random number would be: "
                  << model.random_number_generator()() << "\n";
    }
    return ret;
}

/** Test competencies == 1 with partial table.
 *
 * Values taken from a test without competency table.
 */
int test_two_hosts_with_partial_competency_table_one_only()
{
    int ret = 0;

    Config config;
    config.model_type = "SI";
    config.reproductive_rate = 2;
    config.establishment_probability = 1;
    config.random_seed = 42;
    config.natural_scale = 20;
    config.natural_kernel_type = "deterministic-neighbor";
    config.natural_direction = "E";
    config.use_anthropogenic_kernel = false;
    config.anthro_scale = 0.9;
    config.anthro_kappa = 0;
    config.use_spreadrates = false;
    config.use_quarantine = true;
    config.create_schedules();

    using TestModel = Model<Raster<int>, Raster<double>, Raster<double>::IndexType>;
    TestModel model{config};

    int step = 0;

    Raster<int> total_hosts_1 = {{10, 20, 9}, {0, 0, 0}, {3, 50, 2}};
    Raster<int> total_hosts_2 = {{10, 20, 9}, {14, 15, 0}, {0, 0, 0}};
    Raster<int> infected_1 = {{5, 0, 0}, {0, 0, 0}, {0, 10, 2}};
    Raster<int> infected_2 = {{5, 0, 0}, {0, 5, 0}, {0, 0, 0}};
    Raster<int> total_populations = {{2, 0, 1}, {0, 2, 3}, {2, 2, 5}};
    total_populations += total_hosts_1;
    total_populations += total_hosts_2;
    Raster<int> susceptible_1 = total_hosts_1 - infected_1;
    Raster<int> susceptible_2 = total_hosts_2 - infected_2;
    std::vector<std::vector<int>> suitable_cells =
        find_suitable_cells<Raster<int>::IndexType, Raster<int>>(
            {&total_hosts_1, &total_hosts_2});

    Raster<int> dispersers(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> established_dispersers(infected_1.rows(), infected_1.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    Raster<int> total_exposed_1(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> total_exposed_2(infected_2.rows(), infected_2.cols(), 0);
    Raster<int> died_1(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> died_2(infected_2.rows(), infected_2.cols(), 0);

    Raster<int> empty_integer;
    std::vector<Raster<int>> empty_integers;
    std::vector<Raster<double>> empty_floats;

    unsigned num_mortality_steps = 1;
    std::vector<Raster<int>> mortality_tracker_1(
        num_mortality_steps, Raster<int>(infected_1.rows(), infected_1.cols(), 0));
    std::vector<Raster<int>> mortality_tracker_2(
        num_mortality_steps, Raster<int>(infected_2.rows(), infected_2.cols(), 0));

    std::vector<std::vector<int>> movements = {};

    config.ew_res = 30;
    config.ns_res = 30;
    config.rows = infected_1.rows();
    config.cols = infected_1.cols();

    Raster<double> weather = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    TestModel::StandardSingleHostPool host_pool_1(
        model_type_from_string(config.model_type),
        susceptible_1,
        empty_integers,
        config.latency_period_steps,
        infected_1,
        total_exposed_1,
        empty_integer,
        mortality_tracker_1,
        died_1,
        total_hosts_1,
        model.environment(),
        config.generate_stochasticity,
        config.reproductive_rate,
        config.establishment_stochasticity,
        config.establishment_probability,
        config.rows,
        config.cols,
        suitable_cells);
    TestModel::StandardSingleHostPool host_pool_2(
        model_type_from_string(config.model_type),
        susceptible_2,
        empty_integers,
        config.latency_period_steps,
        infected_2,
        total_exposed_2,
        empty_integer,
        mortality_tracker_2,
        died_2,
        total_hosts_2,
        model.environment(),
        config.generate_stochasticity,
        config.reproductive_rate,
        config.establishment_stochasticity,
        config.establishment_probability,
        config.rows,
        config.cols,
        suitable_cells);
    std::vector<TestModel::StandardSingleHostPool*> host_pools = {
        &host_pool_1, &host_pool_2};
    TestModel::StandardMultiHostPool multi_host_pool(host_pools, config);
    TestModel::StandardPestPool pest_pool{
        dispersers, established_dispersers, outside_dispersers};

    CompetencyTable<TestModel::StandardSingleHostPool> competency_table(
        model.environment());
    competency_table.add_host_competencies({0, 1}, 1);
    competency_table.add_host_competencies({1, 0}, 1);

    multi_host_pool.set_competency_table(competency_table);
    Treatments<TestModel::StandardSingleHostPool, Raster<double>> treatments(
        config.scheduler());
    SpreadRateAction<TestModel::StandardMultiHostPool, int> spread_rate(
        multi_host_pool, config.rows, config.cols, config.ew_res, config.ns_res, 0);
    Raster<int> zeros(infected_1.rows(), infected_1.cols(), 0);
    QuarantineEscapeAction<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, 0);

    model.environment().update_weather_coefficient(weather);

    Raster<int> expected_infected_1 = infected_1;
    Raster<int> expected_infected_2 = infected_2;
    model.run_step(
        step++,
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
    expected_infected_1 += Raster<int>({{0, 4, 0}, {0, 0, 0}, {0, 0, 0}});
    if (infected_1 != expected_infected_1) {
        std::cerr << "test_two_hosts_with_partial_competency_table_one_only (step "
                  << step << ") infected (actual, expected):\n"
                  << infected_1 << "  !=\n"
                  << expected_infected_1 << "\n";
        ++ret;
    }
    expected_infected_2 += Raster<int>({{0, 5, 0}, {0, 0, 0}, {0, 0, 0}});
    if (infected_2 != expected_infected_2) {
        std::cerr << "test_two_hosts_with_partial_competency_table_one_only (step "
                  << step << ") infected (actual, expected):\n"
                  << infected_2 << "  !=\n"
                  << expected_infected_2 << "\n";
        ++ret;
    }
    Raster<int> expected_dispersers = {{16, 0, 0}, {0, 10, 0}, {0, 22, 5}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_two_hosts_with_partial_competency_table_one_only (step "
                  << step << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    size_t expected_outside_dispersers_size = 5;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_two_hosts_with_partial_competency_table_one_only (step "
                  << step << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    model.run_step(
        step++,
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
    expected_infected_1 += Raster<int>({{0, 5, 3}, {0, 0, 0}, {0, 0, 0}});
    if (infected_1 != expected_infected_1) {
        std::cerr << "test_two_hosts_with_partial_competency_table_one_only (step "
                  << step << ") infected (actual, expected):\n"
                  << infected_1 << "  !=\n"
                  << expected_infected_1 << "\n";
        ++ret;
    }
    expected_infected_2 += Raster<int>({{0, 7, 3}, {0, 0, 0}, {0, 0, 0}});
    if (infected_2 != expected_infected_2) {
        std::cerr << "test_two_hosts_with_partial_competency_table_one_only (step "
                  << step << ") infected (actual, expected):\n"
                  << infected_2 << "  !=\n"
                  << expected_infected_2 << "\n";
        ++ret;
    }
    expected_dispersers = {{31, 20, 0}, {0, 13, 0}, {0, 13, 6}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_two_hosts_with_partial_competency_table_one_only (step "
                  << step << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    expected_outside_dispersers_size = 11;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_two_hosts_with_partial_competency_table_one_only (step "
                  << step << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    return ret;
}

/** Test competencies == 1 with complete table.
 *
 * Values taken from a test without competency table.
 */
int test_two_hosts_with_complete_competency_table_one_only()
{
    int ret = 0;

    Config config;
    config.model_type = "SI";
    config.reproductive_rate = 2;
    config.establishment_probability = 1;
    config.random_seed = 42;
    config.natural_scale = 20;
    config.natural_kernel_type = "deterministic-neighbor";
    config.natural_direction = "E";
    config.use_anthropogenic_kernel = false;
    config.anthro_scale = 0.9;
    config.anthro_kappa = 0;
    config.use_spreadrates = false;
    config.use_quarantine = true;
    config.create_schedules();

    using TestModel = Model<Raster<int>, Raster<double>, Raster<double>::IndexType>;
    TestModel model{config};

    int step = 0;

    Raster<int> total_hosts_1 = {{10, 20, 9}, {0, 0, 0}, {3, 50, 2}};
    Raster<int> total_hosts_2 = {{10, 20, 9}, {14, 15, 0}, {0, 0, 0}};
    Raster<int> infected_1 = {{5, 0, 0}, {0, 0, 0}, {0, 10, 2}};
    Raster<int> infected_2 = {{5, 0, 0}, {0, 5, 0}, {0, 0, 0}};
    Raster<int> total_populations = {{2, 0, 1}, {0, 2, 3}, {2, 2, 5}};
    total_populations += total_hosts_1;
    total_populations += total_hosts_2;
    Raster<int> susceptible_1 = total_hosts_1 - infected_1;
    Raster<int> susceptible_2 = total_hosts_2 - infected_2;
    std::vector<std::vector<int>> suitable_cells =
        find_suitable_cells<Raster<int>::IndexType, Raster<int>>(
            {&total_hosts_1, &total_hosts_2});

    Raster<int> dispersers(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> established_dispersers(infected_1.rows(), infected_1.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    Raster<int> total_exposed_1(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> total_exposed_2(infected_2.rows(), infected_2.cols(), 0);
    Raster<int> died_1(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> died_2(infected_2.rows(), infected_2.cols(), 0);

    Raster<int> empty_integer;
    std::vector<Raster<int>> empty_integers;
    std::vector<Raster<double>> empty_floats;

    unsigned num_mortality_steps = 1;
    std::vector<Raster<int>> mortality_tracker_1(
        num_mortality_steps, Raster<int>(infected_1.rows(), infected_1.cols(), 0));
    std::vector<Raster<int>> mortality_tracker_2(
        num_mortality_steps, Raster<int>(infected_2.rows(), infected_2.cols(), 0));

    std::vector<std::vector<int>> movements = {};

    config.ew_res = 30;
    config.ns_res = 30;
    config.rows = infected_1.rows();
    config.cols = infected_1.cols();

    Raster<double> weather = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    TestModel::StandardSingleHostPool host_pool_1(
        model_type_from_string(config.model_type),
        susceptible_1,
        empty_integers,
        config.latency_period_steps,
        infected_1,
        total_exposed_1,
        empty_integer,
        mortality_tracker_1,
        died_1,
        total_hosts_1,
        model.environment(),
        config.generate_stochasticity,
        config.reproductive_rate,
        config.establishment_stochasticity,
        config.establishment_probability,
        config.rows,
        config.cols,
        suitable_cells);
    TestModel::StandardSingleHostPool host_pool_2(
        model_type_from_string(config.model_type),
        susceptible_2,
        empty_integers,
        config.latency_period_steps,
        infected_2,
        total_exposed_2,
        empty_integer,
        mortality_tracker_2,
        died_2,
        total_hosts_2,
        model.environment(),
        config.generate_stochasticity,
        config.reproductive_rate,
        config.establishment_stochasticity,
        config.establishment_probability,
        config.rows,
        config.cols,
        suitable_cells);
    std::vector<TestModel::StandardSingleHostPool*> host_pools = {
        &host_pool_1, &host_pool_2};
    TestModel::StandardMultiHostPool multi_host_pool(host_pools, config);
    TestModel::StandardPestPool pest_pool{
        dispersers, established_dispersers, outside_dispersers};

    std::vector<std::vector<double>> data;
    data.push_back({0, 1, 1});
    data.push_back({1, 0, 1});
    data.push_back({1, 1, 1});
    data.push_back({0, 0, 0});
    config.read_competency_table(data);

    CompetencyTable<TestModel::StandardSingleHostPool> competency_table(
        config, model.environment());

    multi_host_pool.set_competency_table(competency_table);
    Treatments<TestModel::StandardSingleHostPool, Raster<double>> treatments(
        config.scheduler());
    SpreadRateAction<TestModel::StandardMultiHostPool, int> spread_rate(
        multi_host_pool, config.rows, config.cols, config.ew_res, config.ns_res, 0);
    Raster<int> zeros(infected_1.rows(), infected_1.cols(), 0);
    QuarantineEscapeAction<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, 0);

    model.environment().update_weather_coefficient(weather);

    Raster<int> expected_infected_1 = infected_1;
    Raster<int> expected_infected_2 = infected_2;
    model.run_step(
        step++,
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
    expected_infected_1 += Raster<int>({{0, 4, 0}, {0, 0, 0}, {0, 0, 0}});
    if (infected_1 != expected_infected_1) {
        std::cerr << "test_two_hosts_with_complete_competency_table_one_only (step "
                  << step << ") infected (actual, expected):\n"
                  << infected_1 << "  !=\n"
                  << expected_infected_1 << "\n";
        ++ret;
    }
    expected_infected_2 += Raster<int>({{0, 5, 0}, {0, 0, 0}, {0, 0, 0}});
    if (infected_2 != expected_infected_2) {
        std::cerr << "test_two_hosts_with_complete_competency_table_one_only (step "
                  << step << ") infected (actual, expected):\n"
                  << infected_2 << "  !=\n"
                  << expected_infected_2 << "\n";
        ++ret;
    }
    Raster<int> expected_dispersers = {{16, 0, 0}, {0, 10, 0}, {0, 22, 5}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_two_hosts_with_complete_competency_table_one_only (step "
                  << step << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    size_t expected_outside_dispersers_size = 5;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_two_hosts_with_complete_competency_table_one_only (step "
                  << step << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    model.run_step(
        step++,
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
    expected_infected_1 += Raster<int>({{0, 5, 3}, {0, 0, 0}, {0, 0, 0}});
    if (infected_1 != expected_infected_1) {
        std::cerr << "test_two_hosts_with_complete_competency_table_one_only (step "
                  << step << ") infected (actual, expected):\n"
                  << infected_1 << "  !=\n"
                  << expected_infected_1 << "\n";
        ++ret;
    }
    expected_infected_2 += Raster<int>({{0, 7, 3}, {0, 0, 0}, {0, 0, 0}});
    if (infected_2 != expected_infected_2) {
        std::cerr << "test_two_hosts_with_complete_competency_table_one_only (step "
                  << step << ") infected (actual, expected):\n"
                  << infected_2 << "  !=\n"
                  << expected_infected_2 << "\n";
        ++ret;
    }
    expected_dispersers = {{31, 20, 0}, {0, 13, 0}, {0, 13, 6}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_two_hosts_with_complete_competency_table_one_only (step "
                  << step << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    expected_outside_dispersers_size = 11;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_two_hosts_with_complete_competency_table_one_only (step "
                  << step << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    return ret;
}

/** Test competencies (0, 1).
 *
 * Values determined by running the test.
 */
int test_two_hosts_with_table_other_than_one()
{
    int ret = 0;

    Config config;
    config.model_type = "SI";
    config.reproductive_rate = 2;
    config.establishment_probability = 1;
    config.random_seed = 42;
    config.natural_scale = 20;
    config.natural_kernel_type = "deterministic-neighbor";
    config.natural_direction = "E";
    config.use_anthropogenic_kernel = false;
    config.anthro_scale = 0.9;
    config.anthro_kappa = 0;
    config.use_spreadrates = false;
    config.use_quarantine = true;
    config.create_schedules();

    using TestModel = Model<Raster<int>, Raster<double>, Raster<double>::IndexType>;
    TestModel model{config};

    int step = 0;

    Raster<int> total_hosts_1 = {{10, 20, 9}, {0, 0, 0}, {3, 50, 2}};
    Raster<int> total_hosts_2 = {{10, 20, 9}, {14, 15, 0}, {0, 0, 0}};
    Raster<int> infected_1 = {{5, 0, 0}, {0, 0, 0}, {0, 10, 2}};
    Raster<int> infected_2 = {{5, 0, 0}, {0, 5, 0}, {0, 0, 0}};
    Raster<int> total_populations = {{2, 0, 1}, {0, 2, 3}, {2, 2, 5}};
    total_populations += total_hosts_1;
    total_populations += total_hosts_2;
    Raster<int> susceptible_1 = total_hosts_1 - infected_1;
    Raster<int> susceptible_2 = total_hosts_2 - infected_2;
    std::vector<std::vector<int>> suitable_cells =
        find_suitable_cells<Raster<int>::IndexType, Raster<int>>(
            {&total_hosts_1, &total_hosts_2});

    Raster<int> dispersers(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> established_dispersers(infected_1.rows(), infected_1.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    Raster<int> total_exposed_1(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> total_exposed_2(infected_2.rows(), infected_2.cols(), 0);
    Raster<int> died_1(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> died_2(infected_2.rows(), infected_2.cols(), 0);

    Raster<int> empty_integer;
    std::vector<Raster<int>> empty_integers;
    std::vector<Raster<double>> empty_floats;

    unsigned num_mortality_steps = 1;
    std::vector<Raster<int>> mortality_tracker_1(
        num_mortality_steps, Raster<int>(infected_1.rows(), infected_1.cols(), 0));
    std::vector<Raster<int>> mortality_tracker_2(
        num_mortality_steps, Raster<int>(infected_2.rows(), infected_2.cols(), 0));

    std::vector<std::vector<int>> movements = {};

    config.ew_res = 30;
    config.ns_res = 30;
    config.rows = infected_1.rows();
    config.cols = infected_1.cols();

    Raster<double> weather = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    TestModel::StandardSingleHostPool host_pool_1(
        model_type_from_string(config.model_type),
        susceptible_1,
        empty_integers,
        config.latency_period_steps,
        infected_1,
        total_exposed_1,
        empty_integer,
        mortality_tracker_1,
        died_1,
        total_hosts_1,
        model.environment(),
        config.generate_stochasticity,
        config.reproductive_rate,
        config.establishment_stochasticity,
        config.establishment_probability,
        config.rows,
        config.cols,
        suitable_cells);
    TestModel::StandardSingleHostPool host_pool_2(
        model_type_from_string(config.model_type),
        susceptible_2,
        empty_integers,
        config.latency_period_steps,
        infected_2,
        total_exposed_2,
        empty_integer,
        mortality_tracker_2,
        died_2,
        total_hosts_2,
        model.environment(),
        config.generate_stochasticity,
        config.reproductive_rate,
        config.establishment_stochasticity,
        config.establishment_probability,
        config.rows,
        config.cols,
        suitable_cells);
    std::vector<TestModel::StandardSingleHostPool*> host_pools = {
        &host_pool_1, &host_pool_2};
    TestModel::StandardMultiHostPool multi_host_pool(host_pools, config);
    TestModel::StandardPestPool pest_pool{
        dispersers, established_dispersers, outside_dispersers};

    CompetencyTable<TestModel::StandardSingleHostPool> competency_table(
        model.environment());
    competency_table.add_host_competencies({0, 1}, 0.4);
    competency_table.add_host_competencies({1, 0}, 0.6);
    competency_table.add_host_competencies({1, 1}, 0.8);

    multi_host_pool.set_competency_table(competency_table);
    Treatments<TestModel::StandardSingleHostPool, Raster<double>> treatments(
        config.scheduler());
    SpreadRateAction<TestModel::StandardMultiHostPool, int> spread_rate(
        multi_host_pool, config.rows, config.cols, config.ew_res, config.ns_res, 0);
    unsigned quarantine_num_steps =
        get_number_of_scheduled_actions(config.quarantine_schedule());
    Raster<int> zeros(infected_1.rows(), infected_1.cols(), 0);
    QuarantineEscapeAction<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, quarantine_num_steps);

    model.environment().update_weather_coefficient(weather);

    Raster<int> expected_infected_1 = infected_1;
    Raster<int> expected_infected_2 = infected_2;
    model.run_step(
        step++,
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
    expected_infected_1 += Raster<int>({{0, 4, 0}, {0, 0, 0}, {0, 0, 0}});
    if (infected_1 != expected_infected_1) {
        std::cerr << "test_two_hosts_with_table_other_than_one (step " << step
                  << ") infected 1 (actual, expected):\n"
                  << infected_1 << "  !=\n"
                  << expected_infected_1 << "\n";
        ++ret;
    }
    expected_infected_2 += Raster<int>({{0, 4, 0}, {0, 0, 0}, {0, 0, 0}});
    if (infected_2 != expected_infected_2) {
        std::cerr << "test_two_hosts_with_table_other_than_one (step " << step
                  << ") infected 2 (actual, expected):\n"
                  << infected_2 << "  !=\n"
                  << expected_infected_2 << "\n";
        ++ret;
    }
    Raster<int> expected_dispersers = {{14, 0, 0}, {0, 3, 0}, {0, 11, 3}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_two_hosts_with_table_other_than_one (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    size_t expected_outside_dispersers_size = 3;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_two_hosts_with_table_other_than_one (step " << step
                  << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    model.run_step(
        step++,
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
    expected_infected_1 += Raster<int>({{0, 1, 5}, {0, 0, 0}, {0, 0, 0}});
    if (infected_1 != expected_infected_1) {
        std::cerr << "test_two_hosts_with_table_other_than_one (step " << step
                  << ") infected 1 (actual, expected):\n"
                  << infected_1 << "  !=\n"
                  << expected_infected_1 << "\n";
        ++ret;
    }
    expected_infected_2 += Raster<int>({{0, 1, 3}, {0, 0, 0}, {0, 0, 0}});
    if (infected_2 != expected_infected_2) {
        std::cerr << "test_two_hosts_with_table_other_than_one (step " << step
                  << ") infected 2 (actual, expected):\n"
                  << infected_2 << "  !=\n"
                  << expected_infected_2 << "\n";
        ++ret;
    }
    expected_dispersers = {{7, 18, 0}, {0, 3, 0}, {0, 16, 3}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_two_hosts_with_table_other_than_one (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    expected_outside_dispersers_size = 6;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_two_hosts_with_table_other_than_one (step " << step
                  << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    return ret;
}

/** Test two hosts with susceptibilities == 1.
 *
 * Values taken from the test without susceptibilities.
 */
int test_two_hosts_susceptibilities_one()
{
    int ret = 0;

    Config config;
    config.model_type = "SI";
    config.reproductive_rate = 2;
    config.establishment_probability = 1;
    config.random_seed = 42;
    config.natural_scale = 20;
    config.natural_kernel_type = "deterministic-neighbor";
    config.natural_direction = "E";
    config.use_anthropogenic_kernel = false;
    config.anthro_scale = 0.9;
    config.anthro_kappa = 0;
    config.use_spreadrates = false;
    config.use_quarantine = true;
    config.create_schedules();

    using TestModel = Model<Raster<int>, Raster<double>, Raster<double>::IndexType>;
    TestModel model{config};

    int step = 0;

    Raster<int> total_hosts_1 = {{10, 20, 9}, {0, 0, 0}, {3, 50, 2}};
    Raster<int> total_hosts_2 = {{10, 20, 9}, {14, 15, 0}, {0, 0, 0}};
    Raster<int> infected_1 = {{5, 0, 0}, {0, 0, 0}, {0, 10, 2}};
    Raster<int> infected_2 = {{5, 0, 0}, {0, 5, 0}, {0, 0, 0}};
    Raster<int> total_populations = {{2, 0, 1}, {0, 2, 3}, {2, 2, 5}};
    total_populations += total_hosts_1;
    total_populations += total_hosts_2;
    Raster<int> susceptible_1 = total_hosts_1 - infected_1;
    Raster<int> susceptible_2 = total_hosts_2 - infected_2;
    std::vector<std::vector<int>> suitable_cells =
        find_suitable_cells<Raster<int>::IndexType, Raster<int>>(
            {&total_hosts_1, &total_hosts_2});

    Raster<int> dispersers(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> established_dispersers(infected_1.rows(), infected_1.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    Raster<int> total_exposed_1(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> total_exposed_2(infected_2.rows(), infected_2.cols(), 0);
    Raster<int> died_1(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> died_2(infected_2.rows(), infected_2.cols(), 0);

    Raster<int> empty_integer;
    std::vector<Raster<int>> empty_integers;
    std::vector<Raster<double>> empty_floats;

    unsigned num_mortality_steps = 1;
    std::vector<Raster<int>> mortality_tracker_1(
        num_mortality_steps, Raster<int>(infected_1.rows(), infected_1.cols(), 0));
    std::vector<Raster<int>> mortality_tracker_2(
        num_mortality_steps, Raster<int>(infected_2.rows(), infected_2.cols(), 0));

    std::vector<std::vector<int>> movements = {};

    config.ew_res = 30;
    config.ns_res = 30;
    config.rows = infected_1.rows();
    config.cols = infected_1.cols();

    Raster<double> weather = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    TestModel::StandardSingleHostPool host_pool_1(
        model_type_from_string(config.model_type),
        susceptible_1,
        empty_integers,
        config.latency_period_steps,
        infected_1,
        total_exposed_1,
        empty_integer,
        mortality_tracker_1,
        died_1,
        total_hosts_1,
        model.environment(),
        config.generate_stochasticity,
        config.reproductive_rate,
        config.establishment_stochasticity,
        config.establishment_probability,
        config.rows,
        config.cols,
        suitable_cells);
    TestModel::StandardSingleHostPool host_pool_2(
        model_type_from_string(config.model_type),
        susceptible_2,
        empty_integers,
        config.latency_period_steps,
        infected_2,
        total_exposed_2,
        empty_integer,
        mortality_tracker_2,
        died_2,
        total_hosts_2,
        model.environment(),
        config.generate_stochasticity,
        config.reproductive_rate,
        config.establishment_stochasticity,
        config.establishment_probability,
        config.rows,
        config.cols,
        suitable_cells);
    std::vector<TestModel::StandardSingleHostPool*> host_pools = {
        &host_pool_1, &host_pool_2};
    TestModel::StandardMultiHostPool multi_host_pool(host_pools, config);
    TestModel::StandardPestPool pest_pool{
        dispersers, established_dispersers, outside_dispersers};

    PestHostTable<TestModel::StandardSingleHostPool> pest_host_table(
        model.environment());
    pest_host_table.add_host_info(1, 0, 0);
    pest_host_table.add_host_info(1, 0, 0);
    multi_host_pool.set_pest_host_table(pest_host_table);
    Treatments<TestModel::StandardSingleHostPool, Raster<double>> treatments(
        config.scheduler());
    SpreadRateAction<TestModel::StandardMultiHostPool, int> spread_rate(
        multi_host_pool, config.rows, config.cols, config.ew_res, config.ns_res, 0);
    unsigned quarantine_num_steps =
        get_number_of_scheduled_actions(config.quarantine_schedule());
    Raster<int> zeros(infected_1.rows(), infected_1.cols(), 0);
    QuarantineEscapeAction<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, quarantine_num_steps);

    model.environment().update_weather_coefficient(weather);

    Raster<int> expected_infected_1 = infected_1;
    Raster<int> expected_infected_2 = infected_2;
    model.run_step(
        step++,
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
    expected_infected_1 += Raster<int>({{0, 4, 0}, {0, 0, 0}, {0, 0, 0}});
    if (infected_1 != expected_infected_1) {
        std::cerr << "test_two_hosts_susceptibilities_one (step " << step
                  << ") infected (actual, expected):\n"
                  << infected_1 << "  !=\n"
                  << expected_infected_1 << "\n";
        ++ret;
    }
    expected_infected_2 += Raster<int>({{0, 5, 0}, {0, 0, 0}, {0, 0, 0}});
    if (infected_2 != expected_infected_2) {
        std::cerr << "test_two_hosts_susceptibilities_one (step " << step
                  << ") infected (actual, expected):\n"
                  << infected_2 << "  !=\n"
                  << expected_infected_2 << "\n";
        ++ret;
    }
    Raster<int> expected_dispersers = {{16, 0, 0}, {0, 10, 0}, {0, 22, 5}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_two_hosts_susceptibilities_one (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    size_t expected_outside_dispersers_size = 5;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_two_hosts_susceptibilities_one (step " << step
                  << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    model.run_step(
        step++,
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
    expected_infected_1 += Raster<int>({{0, 5, 3}, {0, 0, 0}, {0, 0, 0}});
    if (infected_1 != expected_infected_1) {
        std::cerr << "test_two_hosts_susceptibilities_one (step " << step
                  << ") infected (actual, expected):\n"
                  << infected_1 << "  !=\n"
                  << expected_infected_1 << "\n";
        ++ret;
    }
    expected_infected_2 += Raster<int>({{0, 7, 3}, {0, 0, 0}, {0, 0, 0}});
    if (infected_2 != expected_infected_2) {
        std::cerr << "test_two_hosts_susceptibilities_one (step " << step
                  << ") infected (actual, expected):\n"
                  << infected_2 << "  !=\n"
                  << expected_infected_2 << "\n";
        ++ret;
    }
    expected_dispersers = {{31, 20, 0}, {0, 13, 0}, {0, 13, 6}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_two_hosts_susceptibilities_one (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    expected_outside_dispersers_size = 11;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_two_hosts_susceptibilities_one (step " << step
                  << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    return ret;
}

/** Test two hosts with susceptibilities != 1.
 *
 * Values determined by running the test.
 */
int test_two_hosts_susceptibilities_other_than_one()
{
    int ret = 0;

    Config config;
    config.model_type = "SI";
    config.reproductive_rate = 2;
    config.establishment_probability = 1;
    config.random_seed = 42;
    config.natural_scale = 20;
    config.natural_kernel_type = "deterministic-neighbor";
    config.natural_direction = "E";
    config.use_anthropogenic_kernel = false;
    config.anthro_scale = 0.9;
    config.anthro_kappa = 0;
    config.use_spreadrates = false;
    config.use_quarantine = true;
    config.create_schedules();

    using TestModel = Model<Raster<int>, Raster<double>, Raster<double>::IndexType>;
    TestModel model{config};

    int step = 0;

    Raster<int> total_hosts_1 = {{10, 20, 9}, {0, 0, 0}, {3, 50, 2}};
    Raster<int> total_hosts_2 = {{10, 20, 9}, {14, 15, 0}, {0, 0, 0}};
    Raster<int> infected_1 = {{5, 0, 0}, {0, 0, 0}, {0, 10, 2}};
    Raster<int> infected_2 = {{5, 0, 0}, {0, 5, 0}, {0, 0, 0}};
    Raster<int> total_populations = {{2, 0, 1}, {0, 2, 3}, {2, 2, 5}};
    total_populations += total_hosts_1;
    total_populations += total_hosts_2;
    Raster<int> susceptible_1 = total_hosts_1 - infected_1;
    Raster<int> susceptible_2 = total_hosts_2 - infected_2;
    std::vector<std::vector<int>> suitable_cells =
        find_suitable_cells<Raster<int>::IndexType, Raster<int>>(
            {&total_hosts_1, &total_hosts_2});

    Raster<int> dispersers(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> established_dispersers(infected_1.rows(), infected_1.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;

    Raster<int> total_exposed_1(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> total_exposed_2(infected_2.rows(), infected_2.cols(), 0);
    Raster<int> died_1(infected_1.rows(), infected_1.cols(), 0);
    Raster<int> died_2(infected_2.rows(), infected_2.cols(), 0);

    Raster<int> empty_integer;
    std::vector<Raster<int>> empty_integers;
    std::vector<Raster<double>> empty_floats;

    unsigned num_mortality_steps = 1;
    std::vector<Raster<int>> mortality_tracker_1(
        num_mortality_steps, Raster<int>(infected_1.rows(), infected_1.cols(), 0));
    std::vector<Raster<int>> mortality_tracker_2(
        num_mortality_steps, Raster<int>(infected_2.rows(), infected_2.cols(), 0));

    std::vector<std::vector<int>> movements = {};

    config.ew_res = 30;
    config.ns_res = 30;
    config.rows = infected_1.rows();
    config.cols = infected_1.cols();

    Raster<double> weather = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    TestModel::StandardSingleHostPool host_pool_1(
        model_type_from_string(config.model_type),
        susceptible_1,
        empty_integers,
        config.latency_period_steps,
        infected_1,
        total_exposed_1,
        empty_integer,
        mortality_tracker_1,
        died_1,
        total_hosts_1,
        model.environment(),
        config.generate_stochasticity,
        config.reproductive_rate,
        config.establishment_stochasticity,
        config.establishment_probability,
        config.rows,
        config.cols,
        suitable_cells);
    TestModel::StandardSingleHostPool host_pool_2(
        model_type_from_string(config.model_type),
        susceptible_2,
        empty_integers,
        config.latency_period_steps,
        infected_2,
        total_exposed_2,
        empty_integer,
        mortality_tracker_2,
        died_2,
        total_hosts_2,
        model.environment(),
        config.generate_stochasticity,
        config.reproductive_rate,
        config.establishment_stochasticity,
        config.establishment_probability,
        config.rows,
        config.cols,
        suitable_cells);
    std::vector<TestModel::StandardSingleHostPool*> host_pools = {
        &host_pool_1, &host_pool_2};
    TestModel::StandardMultiHostPool multi_host_pool(host_pools, config);
    TestModel::StandardPestPool pest_pool{
        dispersers, established_dispersers, outside_dispersers};

    PestHostTable<TestModel::StandardSingleHostPool> pest_host_table(
        model.environment());
    pest_host_table.add_host_info(0.8, 0, 0);
    pest_host_table.add_host_info(0.4, 0, 0);
    multi_host_pool.set_pest_host_table(pest_host_table);
    Treatments<TestModel::StandardSingleHostPool, Raster<double>> treatments(
        config.scheduler());
    SpreadRateAction<TestModel::StandardMultiHostPool, int> spread_rate(
        multi_host_pool, config.rows, config.cols, config.ew_res, config.ns_res, 0);
    unsigned quarantine_num_steps =
        get_number_of_scheduled_actions(config.quarantine_schedule());
    Raster<int> zeros(infected_1.rows(), infected_1.cols(), 0);
    QuarantineEscapeAction<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, quarantine_num_steps);

    model.environment().update_weather_coefficient(weather);

    Raster<int> expected_infected_1 = infected_1;
    Raster<int> expected_infected_2 = infected_2;
    model.run_step(
        step++,
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
    expected_infected_1 += Raster<int>({{0, 4, 0}, {0, 0, 0}, {0, 0, 0}});
    if (infected_1 != expected_infected_1) {
        std::cerr << "test_two_hosts_susceptibilities_other_than_one (step " << step
                  << ") infected 1 (actual, expected):\n"
                  << infected_1 << "  !=\n"
                  << expected_infected_1 << "\n";
        ++ret;
    }
    expected_infected_2 += Raster<int>({{0, 3, 0}, {0, 0, 0}, {0, 0, 0}});
    if (infected_2 != expected_infected_2) {
        std::cerr << "test_two_hosts_susceptibilities_other_than_one (step " << step
                  << ") infected 2 (actual, expected):\n"
                  << infected_2 << "  !=\n"
                  << expected_infected_2 << "\n";
        ++ret;
    }
    // First step has the same results as with susceptibility == 1 because
    // number of generated dispersers is influenced by the infected in second step.
    Raster<int> expected_dispersers = {{16, 0, 0}, {0, 10, 0}, {0, 22, 5}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_two_hosts_susceptibilities_other_than_one (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    size_t expected_outside_dispersers_size = 5;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_two_hosts_susceptibilities_other_than_one (step " << step
                  << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    model.run_step(
        step++,
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
    expected_infected_1 += Raster<int>({{0, 4, 1}, {0, 0, 0}, {0, 0, 0}});
    if (infected_1 != expected_infected_1) {
        std::cerr << "test_two_hosts_susceptibilities_other_than_one (step " << step
                  << ") infected 1 (actual, expected):\n"
                  << infected_1 << "  !=\n"
                  << expected_infected_1 << "\n";
        ++ret;
    }
    expected_infected_2 += Raster<int>({{0, 2, 2}, {0, 0, 0}, {0, 0, 0}});
    if (infected_2 != expected_infected_2) {
        std::cerr << "test_two_hosts_susceptibilities_other_than_one (step " << step
                  << ") infected 2 (actual, expected):\n"
                  << infected_2 << "  !=\n"
                  << expected_infected_2 << "\n";
        ++ret;
    }
    expected_dispersers = {{31, 15, 0}, {0, 13, 0}, {0, 17, 1}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_two_hosts_susceptibilities_other_than_one (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    expected_outside_dispersers_size = 6;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_two_hosts_susceptibilities_other_than_one (step " << step
                  << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    return ret;
}

int main()
{
    int ret = 0;

    ret += test_minimal_parameters_one_host();
    ret += test_minimal_parameters_two_hosts();
    ret += test_two_hosts_with_partial_competency_table_one_only();
    ret += test_two_hosts_with_complete_competency_table_one_only();
    ret += test_two_hosts_with_table_other_than_one();
    ret += test_two_hosts_susceptibilities_one();
    ret += test_two_hosts_susceptibilities_other_than_one();
    if (ret) {
        std::cerr << "Test of multi host model: number of errors: " << ret << std::endl;
    }

    return ret;
}

#endif  // POPS_TEST
