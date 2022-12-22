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

int test_update_weather()
{
    int num_errors = 0;
    Environment<Raster<int>, Raster<double>, int> environment;

    Raster<double> step_1{{1, 2}, {0, 1}};
    environment.update_weather_coefficient(step_1);
    if (environment.weather_coefficient_at(0, 0) != 1) {
        std::cerr << "Weather coefficient not set but at-access to it was allowed\n";
        ++num_errors;
    }
    Raster<double> step_2{{1, 2}, {0, 1}};

    return num_errors;
}

int main()
{
    int num_errors = 0;

    num_errors += test_environment_string_values();
    num_errors += test_at_access_rejected();
    num_errors += test_raster_access_rejected();

    std::cout << "Number of errors in environment test: " << num_errors << std::endl;
    return num_errors;
}

#endif  // POPS_TEST
