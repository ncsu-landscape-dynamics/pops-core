/*
 * PoPS model - postprocessing
 *
 * Copyright (C) 2015-2019 by the authors.
 *
 * Authors: Anna Petrasova <akratoc gmail com>
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_POSTPROCESSING_HPP
#define POPS_POSTPROCESSING_HPP

#include <tuple>

namespace pops {

/**
 * Finds bbox of infections and returns
 * north, south, east, west coordinates (as number of rows/cols)
 */
template<typename Raster>
std::tuple<int, int, int, int> infection_boundary(const Raster& raster)
{
    int n = raster.rows() - 1;
    int s = 0;
    int e = 0;
    int w = raster.cols() - 1;
    bool found;
    for (int i = 0; i < raster.rows(); i++) {
        for (int j = 0; j < raster.cols(); j++) {
            auto value = raster(i, j);
            if (value > 0) {
                found = true;
                if (i < n)
                    n = i;
                if (i > s)
                    s = i;
                if (j > e)
                    e = j;
                if (j < w)
                    w = j;
            }
        }
    }
    if (found)
        return std::make_tuple(n, s, e, w);
    else
        return std::make_tuple(-1, -1, -1, -1);
    
}
}
#endif // POPS_POSTPROCESSING_HPP

