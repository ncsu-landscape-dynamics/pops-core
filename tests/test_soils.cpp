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
    std::vector<Raster<int>> rasters{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    Environment<
        Raster<int>,
        Raster<double>,
        Raster<double>::IndexType,
        DefaultSingleGeneratorProvider>
        environment;
    SoilPool<
        Raster<int>,
        Raster<double>,
        Raster<double>::IndexType,
        DefaultSingleGeneratorProvider>
        soils{rasters, environment, false, false, 1};
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
    std::vector<Raster<int>> rasters{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    Environment<
        Raster<int>,
        Raster<double>,
        Raster<double>::IndexType,
        DefaultSingleGeneratorProvider>
        environment;
    SoilPool<
        Raster<int>,
        Raster<double>,
        Raster<double>::IndexType,
        DefaultSingleGeneratorProvider>
        soils{rasters, environment, false, false, 1};
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
 *
 * Values based on the results from the first implementation.
 */
int test_soil_with_model()
{
    int ret = 0;

    Config config;
    config.model_type = "SI";
    config.reproductive_rate = 2;
    config.establishment_probability = 1;
    config.random_seed = 42;
    config.natural_scale = 0.9;
    config.natural_kernel_type = "cauchy";
    config.dispersers_to_soils_percentage = 1.0;
    config.use_anthropogenic_kernel = false;
    config.anthro_scale = 0.9;
    config.anthro_kappa = 0;
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

    config.rows = infected.rows();
    config.cols = infected.cols();

    QuarantineEscapeAction<Raster<int>> quarantine(
        zeros, config.ew_res, config.ns_res, 0);

    Raster<double> weather = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    std::vector<Raster<int>> soil_reservoir(
        2, Raster<int>(infected.rows(), infected.cols(), 0));

    model.environment().update_weather_coefficient(weather);
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
        treatments,
        zeros,
        outside_dispersers,
        quarantine,
        zeros,
        movements,
        Network<int>::null_network(),
        suitable_cells);

    Raster<int> expected_soil_reservoir_0(3, 3, 0);
    Raster<int> expected_soil_reservoir_1 = {{3, 0, 0}, {0, 1, 0}, {0, 0, 0}};
    if (soil_reservoir[0] != expected_soil_reservoir_0) {
        std::cerr << "test_soil_with_model: soil_reservoir[0] (actual, expected):\n"
                  << soil_reservoir[0] << "  !=\n"
                  << expected_soil_reservoir_0 << "\n";
        ++ret;
    }
    if (soil_reservoir[1] != expected_soil_reservoir_1) {
        std::cerr << "test_soil_with_model: soil_reservoir[1] (actual, expected):\n"
                  << soil_reservoir[1] << "  !=\n"
                  << expected_soil_reservoir_1 << "\n";
        ++ret;
    }
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
        treatments,
        zeros,
        outside_dispersers,
        quarantine,
        zeros,
        movements,
        Network<int>::null_network(),
        suitable_cells);
    expected_soil_reservoir_0 = {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    expected_soil_reservoir_1 = {{3, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    if (soil_reservoir[0] != expected_soil_reservoir_0) {
        std::cerr << "test_soil_with_model: soil_reservoir[0] (actual, expected):\n"
                  << soil_reservoir[0] << "  !=\n"
                  << expected_soil_reservoir_0 << "\n";
        ++ret;
    }
    if (soil_reservoir[1] != expected_soil_reservoir_1) {
        std::cerr << "test_soil_with_model: soil_reservoir[1] (actual, expected):\n"
                  << soil_reservoir[1] << "  !=\n"
                  << expected_soil_reservoir_1 << "\n";
        ++ret;
    }
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
