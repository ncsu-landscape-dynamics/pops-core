#ifdef POPS_TEST

/*
 * Simple compilation test for the PoPS Simulation class.
 *
 * Copyright (C) 2018 by the authors.
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

#include "raster.hpp"
#include "simulation.hpp"

#include <map>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>

using std::string;
using std::cout;
using std::cerr;
using std::endl;


int main(int argc, char *argv[])
{
    Raster<int> infected = {{5, 0}, {0, 0}};
    Raster<int> mortality_tracker = {{0, 0}, {0, 0}};
    Raster<int> susceptible = {{5, 6}, {14, 15}};
    Raster<int> total_plants = {{10, 6}, {14, 15}};
    Raster<double> temperature = {{5, 0}, {0, 0}};
    std::vector<std::tuple<int, int>> outside_dispersers;
    Dispersal_kernel dispersal_kernel = CAUCHY;
    double lethal_temperature = -4.5;
    double reproductive_rate = 4.5;
    double short_distance_scale = 18.7;
    Simulation<Raster<int>, Raster<double>> simulation(42, infected);
    simulation.remove(infected, susceptible,
                      temperature, lethal_temperature);
    simulation.generate(infected, 0, reproductive_rate);
    simulation.disperse(susceptible, infected,
                        infected_cohort, total_plants,
                        outside_dispersers, dispersal_kernel, 0,
                        short_distance_scale);
  cout << infected;
  cout << outside_dispersers.size() << endl;
    return 0;
}

#endif  // POPS_TEST
