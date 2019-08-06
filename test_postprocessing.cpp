#ifdef POPS_TEST

/*
 * Simple compilation test for the PoPS postprocessing functions.
 *
 * Copyright (C) 2018-2019 by the authors.
 *
 * Authors: Anna Petrasova <akratoc gmail com>
 *          Vaclav Petras <wenzeslaus gmail com>
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

#include <tuple>
#include <vector>

#include "raster.hpp"
#include "postprocessing.hpp"


using namespace pops;

int test_infected_boundary()
{
    int num_errors = 0;
    Raster<int> infected = {{0, 0, 0, 0, 0},
                            {0, 0, 0, 0, 0},
                            {1, 0, 0, 6, 0},
                            {0, 0, 0, 9, 0},
                            {0, 0, 0, 2, 0}};
    int n, s, e, w;
    std::tie(n, s, e, w) = infection_boundary(infected);

    if (n == 2 && s == 4 && e == 3 && w == 0)
        std::cout << "Infection boundary works (normal case)" << std::endl;

    else {
        std::cout << "Infection boundary does not work  (normal case)" << std::endl;
        std::cout << infected;
        num_errors++;
    }

    Raster<int> infected2 = {{0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 0}};
    std::tie(n, s, e, w) = infection_boundary(infected2);

    if (n == -1) {
        std::cout << "No infection case handled ok" << std::endl;
    }
    else {
        std::cout << "No infection case failed" << std::endl;
        std::cout << infected2;
        num_errors++;
    }
    return num_errors;
}

int test_spread_rate()
{
    int num_errors = 0;
    Raster<int> infected = {{0, 0, 0, 0, 0},
                            {0, 0, 0, 0, 0},
                            {0, 0, 1, 6, 0},
                            {0, 0, 0, 9, 0},
                            {0, 0, 0, 2, 0}};
    Raster<int> infected2 = {{0, 0, 0, 0, 0},
                             {0, 0, 2, 0, 0},
                             {1, 1, 0, 7, 0},
                             {0, 0, 0, 9, 0},
                             {0, 0, 0, 0, 0}};
    bbox_int bbox1 = infection_boundary(infected);
    bbox_int bbox2 = infection_boundary(infected2);
    double n, s, e, w;
    std::tie(n, s, e, w) = spread_rate(bbox1, bbox2, 2, 2, 4);

    if (n == 0.5 && s == -0.5 && e == 0 && w == 1)
        std::cout << "Spread rate works" << std::endl;

    else {
        std::cout << "Spread rate does not work" << std::endl;
        std::cout << n << " " << s << " " << e <<" " << w << std::endl;
        num_errors++;
    }
    return num_errors;
}

int test_spread_rate_average()
{
    int num_errors = 0;
    std::vector<bbox_float> rates;
    rates.push_back(std::make_tuple(0.5, 1, 2, -0.5));
    rates.push_back(std::make_tuple(1, 1, 3, -0.5));
    rates.push_back(std::make_tuple(-0.5, 1, 0, 0.5));
    rates.push_back(std::make_tuple(3, 1, 3, -0.5));
    double n, s, e, w;
    std::tie(n, s, e, w) = average_spread_rate(rates);

    if (n == 1 && s == 1 && e == 2 && w == -0.25)
        std::cout << "Average spread rate works" << std::endl;
    else {
        std::cout << "Average spread rate does not work" << std::endl;
        std::cout << n << " " << s << " " << e <<" " << w << std::endl;
        num_errors++;
    }
    return num_errors;
}

int main()
{
    int num_errors = 0;

    num_errors += test_infected_boundary();
    num_errors += test_spread_rate();
    num_errors += test_spread_rate_average();

    return num_errors;
}

#endif  // POPS_TEST
