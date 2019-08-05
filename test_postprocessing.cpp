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


int main()
{
    int num_errors = 0;

    num_errors += test_infected_boundary();

    return num_errors;
}

#endif  // POPS_TEST
