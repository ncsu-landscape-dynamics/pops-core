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

using namespace pops;

int test_minimal_parameters()
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

    Raster<int> dispersers(infected.rows(), infected.cols());
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

    // std::vector<std::vector<int>> movements = {};
    // Treatments<Raster<int>, Raster<double>> treatments(config.scheduler());
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
    TestModel::StandardMultiHostPool multi_host_pool(host_pools);
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
    model.run_step(
        step++,
        multi_host_pool,
        pest_pool,
        dispersers,
        total_populations,
        treatments,
        empty_floats,
        empty_floats,
        spread_rate,
        quarantine,
        zeros,
        Network<int>::null_network());
    Raster<int> expected_dispersers = {{8, 0, 0}, {0, 8, 0}, {0, 0, 7}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_minimal_parameters (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    size_t expected_outside_dispersers_size = 7;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_minimal_parameters (step " << step
                  << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    model.run_step(
        step++,
        multi_host_pool,
        pest_pool,
        dispersers,
        total_populations,
        treatments,
        empty_floats,
        empty_floats,
        spread_rate,
        quarantine,
        zeros,
        Network<int>::null_network());
    Raster<int> expected_infected = infected;
    if (infected != expected_infected) {
        std::cerr << "test_minimal_parameters (step " << step
                  << ") infected (actual, expected):\n"
                  << infected << "  !=\n"
                  << expected_infected << "\n";
        ++ret;
    }
    expected_dispersers = {{13, 12, 0}, {0, 7, 0}, {0, 0, 1}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_minimal_parameters (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    expected_outside_dispersers_size = 8;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_minimal_parameters (step " << step
                  << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    return ret;
}

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

    Raster<int> dispersers(infected.rows(), infected.cols());
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

    // std::vector<std::vector<int>> movements = {};
    // Treatments<Raster<int>, Raster<double>> treatments(config.scheduler());
    config.ew_res = 30;
    config.ns_res = 30;
    config.rows = infected.rows();
    config.cols = infected.cols();

    Raster<double> weather = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    TestModel::StandardSingleHostPool host_pool_1(
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
    TestModel::StandardSingleHostPool host_pool_2(
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
    std::vector<TestModel::StandardSingleHostPool*> host_pools = {
        &host_pool_1, &host_pool_2};
    TestModel::StandardMultiHostPool multi_host_pool(host_pools);
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
    model.run_step(
        step++,
        multi_host_pool,
        pest_pool,
        dispersers,
        total_populations,
        treatments,
        empty_floats,
        empty_floats,
        spread_rate,
        quarantine,
        zeros,
        Network<int>::null_network());
    Raster<int> expected_dispersers = {{16, 0, 0}, {0, 25, 0}, {0, 0, 6}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    size_t expected_outside_dispersers_size = 6;
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
        dispersers,
        total_populations,
        treatments,
        empty_floats,
        empty_floats,
        spread_rate,
        quarantine,
        zeros,
        Network<int>::null_network());
    Raster<int> expected_infected = infected;
    if (infected != expected_infected) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
                  << ") infected (actual, expected):\n"
                  << infected << "  !=\n"
                  << expected_infected << "\n";
        ++ret;
    }
    expected_dispersers = {{24, 52, 0}, {0, 21, 0}, {0, 0, 6}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    expected_outside_dispersers_size = 12;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
                  << "): outside_dispersers.size (actual, expected): "
                  << outside_dispersers.size()
                  << " != " << expected_outside_dispersers_size << "\n";
        ++ret;
    }
    return ret;
}

int test_minimal_parameters_two_hosts_with_table()
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

    Raster<int> dispersers(infected.rows(), infected.cols());
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

    // std::vector<std::vector<int>> movements = {};
    // Treatments<Raster<int>, Raster<double>> treatments(config.scheduler());
    config.ew_res = 30;
    config.ns_res = 30;
    config.rows = infected.rows();
    config.cols = infected.cols();

    Raster<double> weather = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    TestModel::StandardSingleHostPool host_pool_1(
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
    TestModel::StandardSingleHostPool host_pool_2(
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

    CompetencyTable<TestModel::StandardSingleHostPool, Raster<double>::IndexType>
        competency_table(model.environment());
    competency_table.add_host_competencies({0, 1}, 0.4);
    competency_table.add_host_competencies({1, 0}, 0.1);
    competency_table.add_host_competencies({1, 1}, 0.9);

    std::vector<TestModel::StandardSingleHostPool*> host_pools = {
        &host_pool_1, &host_pool_2};
    TestModel::StandardMultiHostPool multi_host_pool(host_pools);
    TestModel::StandardPestPool pest_pool{
        dispersers, established_dispersers, outside_dispersers};

    multi_host_pool.set_competency_table(competency_table);
    Treatments<TestModel::StandardSingleHostPool, Raster<double>> treatments(
        config.scheduler());
    SpreadRateAction<TestModel::StandardMultiHostPool, int> spread_rate(
        multi_host_pool, config.rows, config.cols, config.ew_res, config.ns_res, 0);
    Raster<int> zeros(infected.rows(), infected.cols(), 0);
    QuarantineEscapeAction<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, 0);

    model.environment().update_weather_coefficient(weather);
    model.run_step(
        step++,
        multi_host_pool,
        pest_pool,
        dispersers,
        total_populations,
        treatments,
        empty_floats,
        empty_floats,
        spread_rate,
        quarantine,
        zeros,
        Network<int>::null_network());
    Raster<int> expected_dispersers = {{16, 0, 0}, {0, 25, 0}, {0, 0, 6}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    size_t expected_outside_dispersers_size = 6;
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
        dispersers,
        total_populations,
        treatments,
        empty_floats,
        empty_floats,
        spread_rate,
        quarantine,
        zeros,
        Network<int>::null_network());
    Raster<int> expected_infected = infected;
    if (infected != expected_infected) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
                  << ") infected (actual, expected):\n"
                  << infected << "  !=\n"
                  << expected_infected << "\n";
        ++ret;
    }
    expected_dispersers = {{24, 52, 0}, {0, 21, 0}, {0, 0, 6}};
    if (dispersers != expected_dispersers) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
                  << ") dispersers (actual, expected):\n"
                  << dispersers << "  !=\n"
                  << expected_dispersers << "\n";
        ++ret;
    }
    expected_outside_dispersers_size = 12;
    if (outside_dispersers.size() != expected_outside_dispersers_size) {
        std::cerr << "test_minimal_parameters_two_hosts (step " << step
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

    ret += test_minimal_parameters();
    ret += test_minimal_parameters_two_hosts();
    ret += test_minimal_parameters_two_hosts_with_table();
    std::cout << "Test of multi host model: number of errors: " << ret << std::endl;

    return ret;
}

#endif  // POPS_TEST
