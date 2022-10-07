#ifdef POPS_TEST

/*
 * Test of PoPS soil borne dispersal.
 *
 * Copyright (C) 2022 by the authors.
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

#include "pops/model.hpp"

using namespace pops;

/**
 * Test basic storage to and retrieval from soils
 */
int test_soils()
{
    int ret = 0;
    std::vector<Raster<int>> rasters{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    Environment<Raster<int>, Raster<double>, Raster<double>::IndexType> environment;
    SoilPool<Raster<int>, Raster<double>, Raster<double>::IndexType> soils{
        rasters, environment, false, false, 1};
    std::default_random_engine generator;

    Raster<double> weather = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
    environment.update_weather_coefficient(weather);

    int initial_dispersers{5};
    soils.dispersers_to(initial_dispersers, 1, 2, generator);
    auto num_dispersers = soils.dispersers_from(1, 2, generator);
    if (num_dispersers != initial_dispersers) {
        std::cerr << "num_dispersers is " << num_dispersers << " (expected all)\n";
        ++ret;
    }
    num_dispersers = soils.dispersers_from(1, 2, generator);
    if (num_dispersers != 0) {
        std::cerr << "num_dispersers is " << num_dispersers << " (expected none)\n";
        ++ret;
    }
    return ret;
}

/**
 * Test influence of weather for soils
 */
int test_soils_weather()
{
    int ret = 0;
    std::vector<Raster<int>> rasters{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    Environment<Raster<int>, Raster<double>, Raster<double>::IndexType> environment;
    SoilPool<Raster<int>, Raster<double>, Raster<double>::IndexType> soils{
        rasters, environment, false, false, 1};
    std::default_random_engine generator;

    Raster<double> weather = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
    double coefficient{0.5};
    weather *= coefficient;
    environment.update_weather_coefficient(weather);

    // The int floor approximation below works only for some even numbers.
    int initial_dispersers{12};
    soils.dispersers_to(initial_dispersers, 1, 2, generator);
    auto num_dispersers = soils.dispersers_from(1, 2, generator);
    auto expected = int(initial_dispersers * coefficient);
    if (num_dispersers != expected) {
        std::cerr << "first num_dispersers is " << num_dispersers << " (expected "
                  << expected << ")\n";
        ++ret;
    }
    num_dispersers = soils.dispersers_from(1, 2, generator);
    expected = int(initial_dispersers * coefficient * coefficient);
    if (num_dispersers != expected) {
        std::cerr << "second num_dispersers is " << num_dispersers << " (expected "
                  << expected << ")\n";
        ++ret;
    }
    num_dispersers = soils.dispersers_from(1, 2, generator);
    expected = int(initial_dispersers * coefficient * coefficient * coefficient);
    if (num_dispersers != expected) {
        std::cerr << "third num_dispersers is " << num_dispersers << " (expected "
                  << expected << ")\n";
        ++ret;
    }
    return ret;
}

/**
 * Test soils runs together with model
 */
int test_soil_with_model()
{
    int ret = 0;

    Config config;
    config.model_type = "SI";
    config.natural_scale = 0.9;
    config.anthro_scale = 0.9;
    config.create_schedules();

    Model<Raster<int>, Raster<double>, Raster<double>::IndexType> model{config};

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

    std::vector<Raster<int>> empty_integer;
    std::vector<Raster<double>> empty_floats;
    Raster<int> zeros(infected.rows(), infected.cols(), 0);

    unsigned num_mortality_steps = 1;
    std::vector<Raster<int>> mortality_tracker(
        num_mortality_steps, Raster<int>(infected.rows(), infected.cols(), 0));

    std::vector<std::vector<int>> movements = {};
    Treatments<Raster<int>, Raster<double>> treatments(config.scheduler());
    config.ew_res = 30;
    config.ns_res = 30;
    unsigned rate_num_steps =
        get_number_of_scheduled_actions(config.spread_rate_schedule());
    SpreadRate<Raster<int>> spread_rate(
        infected, config.ew_res, config.ns_res, rate_num_steps, suitable_cells);
    QuarantineEscape<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, 0, suitable_cells);

    Raster<double> weather = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    std::vector<Raster<int>> soil_reservoir(
        1, Raster<int>(infected.rows(), infected.cols(), 0));

    model.activate_soils(soil_reservoir);
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
        weather,
        treatments,
        zeros,
        outside_dispersers,
        spread_rate,
        quarantine,
        zeros,
        movements,
        Network<int>::null_network(),
        suitable_cells);
    return ret;
}

int main()
{
    int ret = 0;

    ret += test_soils();
    ret += test_soils_weather();
    ret += test_soil_with_model();

    return ret;
}

#endif  // POPS_TEST
