#ifdef POPS_TEST

/*
 * Tst for the PoPS Environment class.
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

#include <vector>
#include <tuple>

#include <pops/environment.hpp>
#include <pops/raster.hpp>

using namespace pops;

template<typename Exception, typename Operation>
bool throws_exception(Operation operation)
{
    try {
        operation();
    }
    catch (const Exception&) {
        return true;
    }
    return false;
}

int test_environment_string_values()
{
    int num_errors = 0;
    weather_type_from_string("deterministic");
    weather_type_from_string("Deterministic");
    weather_type_from_string("probabilistic");
    weather_type_from_string("Probabilistic");
    weather_type_from_string("");
    weather_type_from_string("none");
    weather_type_from_string("None");
    weather_type_from_string("NONE");
    for (const auto& value : {"does-not-exist", "PROBABILISTIC", "prob"}) {
        bool thrown = throws_exception<std::invalid_argument>(
            [&value] { weather_type_from_string(value); });
        if (!thrown) {
            std::cerr << "Value '" << value
                      << "' wrongly accepted by weather_type_from_string()\n";
            ++num_errors;
        }
    }
    return num_errors;
}

int test_at_access_rejected()
{
    Environment<Raster<int>, Raster<double>, int> environment;
    try {
        auto a = environment.weather_coefficient_at(0, 0);
        std::cout << a;
    }
    catch (const std::logic_error&) {
        return 0;
    }
    std::cerr << "Weather coefficient not set but at-access to it was allowed\n";
    return 1;
}

int test_raster_access_rejected()
{
    Environment<Raster<int>, Raster<double>, int> environment;
    try {
        auto& a = environment.weather_coefficient();
        std::cout << a;
    }
    catch (const std::logic_error&) {
        return 0;
    }
    std::cerr << "Weather coefficient not set but at-access to it was allowed\n";
    return 1;
}

int test_update_deterministic_weather()
{
    int num_errors = 0;
    Environment<Raster<int>, Raster<double>, int> environment;

    Raster<double> step_1{{1, 2}, {0, 1}, {8, 7}};
    environment.update_weather_coefficient(step_1);
    for (int row = 0; row < step_1.rows(); ++row) {
        for (int col = 0; col < step_1.cols(); ++col) {
            if (environment.weather_coefficient_at(row, col) != step_1(row, col)) {
                std::cerr << "Deterministic weather coefficient step 1 at (" << row
                          << "," << col << ") equals to "
                          << environment.weather_coefficient_at(row, col)
                          << " but it should be " << step_1(row, col) << "\n";
                ++num_errors;
            }
        }
    }
    Raster<double> step_2{{0, 5}, {1, 2}, {4, 3}};
    environment.update_weather_coefficient(step_2);
    for (int row = 0; row < step_2.rows(); ++row) {
        for (int col = 0; col < step_2.cols(); ++col) {
            if (environment.weather_coefficient_at(row, col) != step_2(row, col)) {
                std::cerr << "Deterministic weather coefficient step 2 at (" << row
                          << "," << col << ") equals to "
                          << environment.weather_coefficient_at(row, col)
                          << " but it should be " << step_2(row, col) << "\n";
                ++num_errors;
            }
        }
    }

    return num_errors;
}

int test_update_probabilistic_weather()
{
    int num_errors = 0;
    Environment<Raster<int>, Raster<double>, int> environment;
    unsigned seed = 1;
    std::default_random_engine generator(seed);

    Raster<double> step_1{{0.1, 0.2}, {0.0, 0.3}, {0.8, 0.7}};
    Raster<double> stddev_step_1(step_1.rows(), step_1.cols(), 0);
    environment.update_weather_from_distribution(step_1, stddev_step_1, generator);
    for (int row = 0; row < step_1.rows(); ++row) {
        for (int col = 0; col < step_1.cols(); ++col) {
            if (environment.weather_coefficient_at(row, col) != step_1(row, col)) {
                std::cerr << "Probabilistic weather coefficient step 1 at (" << row
                          << "," << col << ") equals to "
                          << environment.weather_coefficient_at(row, col)
                          << " but it should be " << step_1(row, col) << "\n";
                ++num_errors;
            }
        }
    }
    // Semi-lenient test which works with seeds 0, 1, 5, 6, 7.
    double tolerance = 0.6;
    Raster<double> step_2{{0, 0.5}, {0.1, 0.2}, {0.4, 0.3}};
    Raster<double> stddev_step_2{{0.5, 5.8}, {0.1, 2}, {4, 3}};
    environment.update_weather_from_distribution(step_2, stddev_step_2, generator);
    for (int row = 0; row < step_2.rows(); ++row) {
        for (int col = 0; col < step_2.cols(); ++col) {
            auto difference =
                environment.weather_coefficient_at(row, col) - step_2(row, col);
            if (difference > tolerance) {
                std::cerr << "Probabilistic weather coefficient step 2 at (" << row
                          << "," << col << ") equals to "
                          << environment.weather_coefficient_at(row, col)
                          << " but it should closer to " << step_2(row, col)
                          << " (difference " << difference
                          << " is greater than tolerance " << tolerance << ")\n";
                ++num_errors;
            }
        }
    }
    return num_errors;
}

/** The other tests uses more rows than columns, so this one just checks the case with
 * more columns than rows.
 */
int test_update_probabilistic_weather_dimensions()
{
    int num_errors = 0;
    Environment<Raster<int>, Raster<double>, int> environment;
    std::default_random_engine generator(1);

    Raster<double> step_1{{0.1, 0.2, 0.0, 0.3, 0.8, 0.7}};
    Raster<double> stddev_step_1(step_1.rows(), step_1.cols(), 0);
    environment.update_weather_from_distribution(step_1, stddev_step_1, generator);
    for (int row = 0; row < step_1.rows(); ++row) {
        for (int col = 0; col < step_1.cols(); ++col) {
            if (environment.weather_coefficient_at(row, col) != step_1(row, col)) {
                std::cerr << "Probabilistic weather coefficient with many columns at ("
                          << row << "," << col << ") equals to "
                          << environment.weather_coefficient_at(row, col)
                          << " but it should be " << step_1(row, col) << "\n";
                ++num_errors;
            }
        }
    }
    return num_errors;
}

int test_update_probabilistic_weather_range()
{
    int num_errors = 0;
    Environment<Raster<int>, Raster<double>, int> environment;
    std::default_random_engine generator(1);
    int num_tests = 10;

    Raster<double> mean{{0, 0.5}, {1, 0.2}, {0.4, 0.3}};
    Raster<double> stddev{{0.5, 5.8}, {10, 20}, {4, 30}};
    for (int i = 0; i < num_tests; ++i) {
        environment.update_weather_from_distribution(mean, stddev, generator);
        for (int row = 0; row < mean.rows(); ++row) {
            for (int col = 0; col < mean.cols(); ++col) {
                auto coeff = environment.weather_coefficient_at(row, col);
                if (coeff < 0) {
                    std::cerr << "Probabilistic weather coefficient lower than zero: "
                              << coeff << " (at " << row << "," << col << ")\n";
                    ++num_errors;
                }
                else if (coeff > 1) {
                    std::cerr << "Probabilistic weather coefficient greater than one: "
                              << coeff << " (at " << row << "," << col << ")\n";
                    ++num_errors;
                }
            }
        }
    }
    return num_errors;
}

int test_update_probabilistic_weather_mean_in_range()
{
    int num_errors = 0;
    Environment<Raster<int>, Raster<double>, int> environment;
    unsigned seed = 1;
    std::default_random_engine generator(seed);

    Raster<double> step_1{{0.1, 0.2}, {10.0, 0.3}, {0.8, 0.7}};
    Raster<double> stddev_step_1(step_1.rows(), step_1.cols(), 0);

    bool thrown = throws_exception<std::invalid_argument>(
        [&environment, &step_1, &stddev_step_1, &generator] {
            environment.update_weather_from_distribution(
                step_1, stddev_step_1, generator);
        });
    if (!thrown) {
        std::cerr << "An out-of-range value wrongly accepted "
                  << "by update_weather_from_distribution()\n";
        ++num_errors;
    }
    return num_errors;
}

int main()
{
    int num_errors = 0;

    num_errors += test_environment_string_values();
    num_errors += test_at_access_rejected();
    num_errors += test_raster_access_rejected();
    num_errors += test_update_deterministic_weather();
    num_errors += test_update_probabilistic_weather();
    num_errors += test_update_probabilistic_weather_dimensions();
    num_errors += test_update_probabilistic_weather_range();

    std::cout << "Number of errors in environment test: " << num_errors << std::endl;
    return num_errors;
}

#endif  // POPS_TEST
