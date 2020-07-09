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

#include "model.hpp"

using namespace pops;
using std::cout;

int test_with_reduced_stochasticity()
{
    Raster<int> infected = {{5, 0}, {0, 0}};
    Raster<int> susceptible = {{10, 20}, {14, 15}};
    Raster<int> total_hosts = susceptible;
    Raster<int> zeros(infected.rows(), infected.cols(), 0);

    Raster<int> expected_mortality_tracker = {{0, 10}, {0, 0}};
    auto expected_infected = expected_mortality_tracker + infected;

    Raster<int> dispersers(infected.rows(), infected.cols());
    std::vector<std::tuple<int, int>> outside_dispersers;
    Config config;
    config.weather = false;
    config.reproductive_rate = 2;
    config.generate_stochasticity = false;
    config.establishment_stochasticity = false;
    // We want everything to establish.
    config.establishment_probability = 1;
    config.natural_kernel_type = "deterministic_neighbor";
    config.natural_direction = "E";
    config.use_anthropogenic_kernel = false;
    config.random_seed = 42;
    config.rows = infected.rows();
    config.cols = infected.cols();
    config.model_type = "SI";
    config.latency_period_steps = 0;
    config.use_lethal_temperature = false;

    config.set_date_start(2020, 1, 1);
    config.set_date_end(2021, 12, 31);
    config.set_step_unit(StepUnit::Month);
    config.set_step_num_units(1);
    config.create_schedules();

    int weather_step = 0;
    unsigned num_mortality_years = config.num_mortality_years();
    std::cerr << "num_mortality_years: " << num_mortality_years << "\n";
    std::vector<Raster<int>> mortality_tracker(
        num_mortality_years, Raster<int>(infected.rows(), infected.cols(), 0));

    //    int exposed_size = 0;
    //    if (config.latency_period_steps)
    //        exposed_size = config.latency_period_steps + 1;
    //    std::vector<Raster<int>> exposed(
    //                exposed_size,
    //                Raster<int>(infected.rows(), infected.cols(), 0));
    Raster<int> died(infected.rows(), infected.cols(), 0);
    std::vector<Raster<int>> empty_integer;
    std::vector<Raster<double>> empty_float;
    Treatments<Raster<int>, Raster<double>> treatments(config.scheduler());
    config.use_treatments = false;
    config.ew_res = 1;
    config.ns_res = 1;
    unsigned rate_num_years =
        get_number_of_scheduled_actions(config.spread_rate_schedule());
    SpreadRate<Raster<int>> spread_rate(
        infected, config.ew_res, config.ns_res, rate_num_years);

    auto expected_dispersers = config.reproductive_rate * infected;

    int step = 0;

    Model<Raster<int>, Raster<double>, Raster<double>::IndexType> model(config);
    model.run_step(
        step++,
        weather_step,
        infected,
        susceptible,
        total_hosts,
        dispersers,
        empty_integer,
        mortality_tracker,
        died,
        empty_float,
        empty_float,
        treatments,
        zeros,
        outside_dispersers,
        spread_rate);
    if (dispersers != expected_dispersers) {
        cout << "reduced_stochasticity: dispersers (actual, expected):\n"
             << dispersers << "  !=\n"
             << expected_dispersers << "\n";
        return 1;
    }
    if (!outside_dispersers.empty()) {
        cout << "reduced_stochasticity: There are outside_dispersers ("
             << outside_dispersers.size() << ") but there should be none\n";
        return 1;
    }
    if (infected != expected_infected) {
        cout << "reduced_stochasticity: infected (actual, expected):\n"
             << infected << "  !=\n"
             << expected_infected << "\n";
        return 1;
    }
    if (mortality_tracker[0] != expected_mortality_tracker) {
        cout << "reduced_stochasticity: mortality tracker (actual, expected):\n"
             << mortality_tracker[0] << "  !=\n"
             << expected_mortality_tracker << "\n";
        return 1;
    }
    return 0;
}

int main()
{
    int ret = 0;

    ret += test_with_reduced_stochasticity();

    return ret;
}

#endif  // POPS_TEST
