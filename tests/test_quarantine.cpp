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


#include "raster.hpp"
#include "quarantine.hpp"

using namespace pops;

int test_quarantine()
{
    int err = 0;
    Raster<int> areas = {
        {0, 1, 1, 0, 0},
        {0, 0, 1, 4, 0},
        {0, 0, 0, 4, 0},
        {0, 3, 0, 4, 4},
        {0, 0, 0, 4, 0}};

    Raster<int> infectionOutside = {
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0}};
    
    Raster<int> infectionInside = {
        {0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 0}};

    QuarantineEscape<Raster<int>> quarantine(areas, 10, 10);
    EscapeDistDir res_out = quarantine.infection_inside_quarantine(infectionOutside, areas);
    EscapeDistDir res_in = quarantine.infection_inside_quarantine(infectionInside, areas);
    std::vector<DistDir> distdir;
    bool escaped;
    double dist;
    Direction d;
    
    std::tie(escaped, distdir) = res_out;
    std::cout << escaped << std::endl;
    for (unsigned i = 0; i < distdir.size(); i++) {
        std::tie(dist, d) = distdir[i];
        std::cout << "Min dist: " << dist << ", direction: " << d << std::endl;
    }
    
    std::tie(escaped, distdir) = res_in;
    std::cout << escaped << std::endl;
    for (unsigned i = 0; i < distdir.size(); i++) {
        std::tie(dist, d) = distdir[i];
        std::cout << "Min dist: " << dist << ", direction: " << d << std::endl;
    }
    
//    double n, s, e, w;
//    std::tie(n, s, e, w) = spread_rate.yearly_rate(0);
//    if (!(n == 0 && s == 0 && e == 10 && w == 10)) {
//        std::cout << "spread rate for year 1 fails" << std::endl;
//        err++;
//    }

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
