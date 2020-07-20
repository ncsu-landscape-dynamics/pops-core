#ifdef POPS_TEST

/*
 * Simple compilation test for the PoPS quarantine class.
 *
 * Copyright (C) 2020 by the authors.
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

#include <pops/raster.hpp>
#include <pops/quarantine.hpp>

using namespace pops;

int test_quarantine()
{
    int err = 0;
    Raster<int> areas = {
        {0, 1, 1, 0, 0},
        {0, 0, 1, 4, 0},
        {0, 0, 4, 4, 4},
        {0, 3, 4, 4, 4},
        {0, 0, 0, 4, 0}};

    // infected rasters for 1st run and 3 years
    Raster<int> infection11 = {
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0}};

    Raster<int> infection12 = {
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 0}};

    Raster<int> infection13 = {
        {0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 1, 1, 0},
        {0, 0, 0, 0, 0}};

    // infected rasters for 2nd run and 3 years
    Raster<int> infection21 = {
        {0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0}};

    Raster<int> infection22 = {
        {0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0}};

    Raster<int> infection23 = {
        {0, 0, 0, 1, 0},
        {0, 0, 1, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 1, 1, 0},
        {0, 0, 0, 0, 0}};

    std::vector<QuarantineEscape<Raster<int>>> runs(
        2, QuarantineEscape<Raster<int>>(areas, 10, 10, 3));
    runs[0].infection_escape_quarantine(infection11, areas, 0);
    runs[0].infection_escape_quarantine(infection12, areas, 1);
    runs[0].infection_escape_quarantine(infection13, areas, 2);
    runs[1].infection_escape_quarantine(infection21, areas, 0);
    runs[1].infection_escape_quarantine(infection22, areas, 1);
    runs[1].infection_escape_quarantine(infection23, areas, 2);
    double p1 = quarantine_escape_probability(runs, 0);
    double p2 = quarantine_escape_probability(runs, 1);
    double p3 = quarantine_escape_probability(runs, 2);
    if (!(p1 == 0 && p2 == 0 && p3 == 0.5)) {
        std::cout << "Probability of escape fails" << std::endl;
        std::cout << p1 << " " << p2 << " " << p3 << std::endl;
        err++;
    }

    std::vector<double> d1 = distance_to_quarantine(runs, 0);
    std::vector<double> d2 = distance_to_quarantine(runs, 1);
    std::vector<double> d3 = distance_to_quarantine(runs, 2);
    if (!(d1.at(0) == 10 && d1.at(1) == 0)) {
        std::cout << "Distance to quarantine boundary fails" << std::endl;
        std::cout << d1.at(0) << " " << d1.at(1) << std::endl;
        err++;
    }
    if (!(d2.at(0) == 10 && d2.at(1) == 0)) {
        std::cout << "Distance to quarantine boundary fails" << std::endl;
        std::cout << d2.at(0) << " " << d2.at(1) << std::endl;
        err++;
    }
    if (!(d3.at(0) == 0 && std::isnan(d3.at(1)))) {
        std::cout << "Distance to quarantine boundary fails" << std::endl;
        std::cout << d3.at(0) << " " << d3.at(1) << std::endl;
        err++;
    }

    // just test calling the function
    std::string output = write_quarantine_escape(runs, 3, 2000);
    if (output.empty())
        err++;

    return err;
}

int main()
{
    int num_errors = 0;

    num_errors += test_quarantine();
    std::cout << "Quarantine number of errors: " << num_errors << std::endl;
    return num_errors;
}

#endif  // POPS_TEST
