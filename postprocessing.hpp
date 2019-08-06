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
#include <vector>

namespace pops {

typedef std::tuple<int, int, int, int> bbox_int;
typedef std::tuple<double, double, double, double> bbox_float;

/**
 * Finds bbox of infections and returns
 * north, south, east, west coordinates (as number of rows/cols)
 */
template<typename Raster>
bbox_int infection_boundary(const Raster& raster)
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

/**
 * Computes spread rate in n, s, e, w directions
 * based on bounding boxes of disease (bbox2 - bbox1).
 * Unit is distance (map units) per year.
 */
bbox_float spread_rate(bbox_int bbox1, bbox_int bbox2, double ew_res, double ns_res, int years = 1)
{
    int n1, n2, s1, s2, e1, e2, w1, w2;
    std::tie(n1, s1, e1, w1) = bbox1;
    std::tie(n2, s2, e2, w2) = bbox2;
    double n_rate = ((n1 - n2) * ns_res) / years;
    double s_rate = ((s2 - s1) * ns_res) / years;
    double e_rate = ((e2 - e1) * ew_res) / years;
    double w_rate = ((w1 - w2) * ew_res) / years;

    return std::make_tuple(n_rate, s_rate, e_rate, w_rate);
}

/**
 * Computes average spread rate in n, s, e, w directions
 * from vector of spread rates.
 */
bbox_float average_spread_rate(std::vector<bbox_float> bboxes)
{
    int size = bboxes.size();
    double n, s, e, w;
    double avg_n = 0, avg_s = 0, avg_e = 0, avg_w = 0;
    for (int i = 0; i < size; i++) {
        std::tie(n, s, e, w) = bboxes[i];
        avg_n += n;
        avg_s += s;
        avg_e += e;
        avg_w += w;
    }
    avg_n /= size;
    avg_s /= size;
    avg_e /= size;
    avg_w /= size;

    return std::make_tuple(avg_n, avg_s, avg_e, avg_w);
}

}
#endif // POPS_POSTPROCESSING_HPP

