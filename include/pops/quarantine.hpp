/*
 * PoPS model - Quarantine escape computation
 *
 * Copyright (C) 2020 by the authors.
 *
 * Authors: Anna Petrasova <akratoc ncsu edu>
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_QUARANTINE_HPP
#define POPS_QUARANTINE_HPP

#include <tuple>
#include <map>
#include <vector>
#include <cmath>
#include <limits>
#include <type_traits>

namespace pops {

/*! Quarantine direction
 */
enum class Direction
{
    N = 0,  //!< North
    E = 90,  //!< NEast
    S = 180,  //!< South
    W = 270  //!< West
};
std::ostream& operator<<(std::ostream& os, const Direction& obj)
{
    os << static_cast<std::underlying_type<Direction>::type>(obj);
    return os;
}

typedef std::tuple<int, int, int, int> BBoxInt;
typedef std::tuple<double, double, double, double> BBoxFloat;
typedef std::tuple<bool, bool, bool, bool> BBoxBool;
typedef std::tuple<double, Direction> DistDir;
typedef std::tuple<bool, std::vector<DistDir>> EscapeDistDir;

template<typename IntegerRaster, typename RasterIndex = int>
class QuarantineEscape
{
private:
    RasterIndex width_;
    RasterIndex height_;
    // the west-east resolution of the pixel
    double west_east_resolution;
    // the north-south resolution of the pixel
    double north_south_resolution;
    std::vector<BBoxInt> boundaries;
    std::map<int, int> boundary_id_idx_map;
    IntegerRaster quarantine_areas_;

    void quarantine_boundary(const IntegerRaster& quarantine_areas)
    {
        int n, s, e, w;
        int idx = 0;
        for (int i = 0; i < height_; i++) {
            for (int j = 0; j < width_; j++) {
                auto value = quarantine_areas(i, j);
                if (value > 0) {
                    auto search = boundary_id_idx_map.find(value);
                    int bidx;
                    if (search == boundary_id_idx_map.end()) {
                        boundary_id_idx_map.insert(std::make_pair(value, idx));
                        boundaries.push_back(
                            std::make_tuple(height_ - 1, 0, 0, width_ - 1));
                        bidx = idx;
                        ++idx;
                    }
                    else
                        bidx = search->second;
                    std::tie(n, s, e, w) = boundaries.at(bidx);
                    if (i < n)
                        n = i;
                    if (i > s)
                        s = i;
                    if (j > e)
                        e = j;
                    if (j < w)
                        w = j;
                    boundaries.at(bidx) = std::make_tuple(n, s, e, w);
                }
            }
        }
    }
    std::tuple<double, Direction>
    closest_direction(RasterIndex i, RasterIndex j, const BBoxInt boundary)
    {
        int n, s, e, w;
        int mindist = std::numeric_limits<int>::max();
        std::tie(n, s, e, w) = boundary;
        std::cout << n << " " << s << " " << e << " " << w << std::endl;
        std::cout << i << " " << j << std::endl;
        std::tuple<double, Direction> closest;
        if ((i - n) * north_south_resolution < mindist) {
            mindist = (i - n) * north_south_resolution;
            closest = std::make_tuple(mindist, Direction::N);
        }
        if ((s - i) * north_south_resolution < mindist) {
            mindist = (s - i) * north_south_resolution;
            closest = std::make_tuple(mindist, Direction::S);
        }
        if ((e - j) * west_east_resolution < mindist) {
            mindist = (e - j) * west_east_resolution;
            closest = std::make_tuple(mindist, Direction::E);
        }
        if ((j - w) * west_east_resolution < mindist) {
            mindist = (j - w) * west_east_resolution;
            closest = std::make_tuple(mindist, Direction::W);
        }
        std::cout << mindist << std::endl;
        return closest;
    }

public:
    QuarantineEscape(
        const IntegerRaster& quarantine_areas, double ew_res, double ns_res)
        : width_(quarantine_areas.cols()),
          height_(quarantine_areas.rows()),
          west_east_resolution(ew_res),
          north_south_resolution(ns_res)
    {
        quarantine_boundary(quarantine_areas);
        double n, s, e, w;
        for (unsigned i = 0; i < boundaries.size(); i++) {
            std::tie(n, s, e, w) = boundaries.at(i);
        }
    }

    QuarantineEscape() = delete;

    EscapeDistDir infection_inside_quarantine(
        const IntegerRaster& infected, const IntegerRaster& quarantine_areas)
    {
        std::vector<DistDir> min_dist_dir;
        for (unsigned i = 0; i < boundaries.size(); i++) {
            min_dist_dir.push_back(
                std::make_tuple(std::numeric_limits<int>::max(), Direction::N));
        }
        for (int i = 0; i < height_; i++) {
            for (int j = 0; j < width_; j++) {
                if (!infected(i, j))
                    continue;
                int area = quarantine_areas(i, j);
                if (area == 0) {
                    return std::make_tuple(false, min_dist_dir);
                }
                double dist;
                Direction dir;
                int bindex = boundary_id_idx_map[area];
                std::tie(dist, dir) = closest_direction(i, j, boundaries.at(bindex));
                if (dist < std::get<0>(min_dist_dir.at(bindex))) {
                    min_dist_dir.at(bindex) = std::make_tuple(dist, dir);
                }
            }
        }
        return std::make_tuple(true, min_dist_dir);
    }
};
}  // namespace pops
#endif  // POPS_QUARANTINE_HPP
