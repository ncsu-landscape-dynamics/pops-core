#ifdef POPS_TEST

/*
 * Tests if random number generators appropriately approximates the power law
 * distribution and the hyperbolic secant distribution
 *
 * Copyright (C) 2023 by the authors.
 *
 * Authors: Vaclav Petras (wenzeslaus gmail com)
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

#include <iostream>

#include <pops/normal_distribution_with_uniform_fallback.hpp>

using std::string;
using std::cout;

using namespace pops;

/**
 * Check that the distribution produces values in a specified range.
 */
int test_normal_with_uniform_fallback()
{
    int num_errors = 0;
    unsigned seed = 1;
    std::default_random_engine generator(seed);
    double low = 11;
    double high = 12;
    NormalDistributionWithUniformFallback<double> distribution(11.8, 2, low, high);

    for (int i = 0; i < 100000000; i++) {
        double x = distribution(generator);
        if (x < low || x > high) {
            std::cerr << x << " is out-of-range [" << low << ", " << high << "]\n";
            ++num_errors;
        }
    }
    return num_errors;
}
int main()
{
    int num_errors = 0;

    num_errors += test_normal_with_uniform_fallback();

    std::cout << "Distributions test number of errors: " << num_errors << "\n";
    return num_errors;
}
#endif  // POPS_TEST
