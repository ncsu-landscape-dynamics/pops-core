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
#include <sstream>
#include <iomanip>

namespace pops {

/*! Quarantine direction
 */
enum class Direction
{
    N = 0,  //!< North
    E = 90,  //!< East
    S = 180,  //!< South
    W = 270,  //!< West
    None  //!< Escaped
};
std::ostream& operator<<(std::ostream& os, const Direction& obj)
{
    os << static_cast<std::underlying_type<Direction>::type>(obj);
    return os;
}

typedef std::tuple<int, int, int, int> BBoxInt;
typedef std::tuple<double, Direction> DistDir;
typedef std::tuple<bool, DistDir> EscapeDistDir;
typedef std::vector<EscapeDistDir> EscapeDistDirs;

/**
 * Class storing and computing quarantine escap metrics for one simulation.
 */
template<typename IntegerRaster, typename RasterIndex = int>
class QuarantineEscape
{
private:
    RasterIndex width_;
    RasterIndex height_;
    // the west-east resolution of the pixel
    double west_east_resolution_;
    // the north-south resolution of the pixel
    double north_south_resolution_;
    unsigned num_years;
    std::vector<BBoxInt> boundaries;
    // mapping between quarantine areas is from map and index
    std::map<int, int> boundary_id_idx_map;
    std::vector<EscapeDistDir> escape_dist_dirs;

    /**
     * Computes bbox of each quarantine area.
     * Different quarantine areas are represented by different integers.
     * 0 in the raster means no quarantine area.
     */
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
    /**
     * Computes minimum distance (in map units) and the associated direction
     * to quarantine area boundary.
     * @param i infected cell row
     * @param j infected cell col
     * @param boundary quarantine area boundary
     */

    DistDir
    closest_direction(RasterIndex i, RasterIndex j, const BBoxInt boundary) const
    {
        int n, s, e, w;
        int mindist = std::numeric_limits<int>::max();
        std::tie(n, s, e, w) = boundary;
        DistDir closest;
        if ((i - n) * north_south_resolution_ < mindist) {
            mindist = (i - n) * north_south_resolution_;
            closest = std::make_tuple(mindist, Direction::N);
        }
        if ((s - i) * north_south_resolution_ < mindist) {
            mindist = (s - i) * north_south_resolution_;
            closest = std::make_tuple(mindist, Direction::S);
        }
        if ((e - j) * west_east_resolution_ < mindist) {
            mindist = (e - j) * west_east_resolution_;
            closest = std::make_tuple(mindist, Direction::E);
        }
        if ((j - w) * west_east_resolution_ < mindist) {
            mindist = (j - w) * west_east_resolution_;
            closest = std::make_tuple(mindist, Direction::W);
        }
        return closest;
    }

public:
    QuarantineEscape(
        const IntegerRaster& quarantine_areas,
        double ew_res,
        double ns_res,
        unsigned num_years)
        : width_(quarantine_areas.cols()),
          height_(quarantine_areas.rows()),
          west_east_resolution_(ew_res),
          north_south_resolution_(ns_res),
          num_years(num_years),
          escape_dist_dirs(
              num_years,
              std::make_tuple(
                  false,
                  std::make_tuple(std::numeric_limits<double>::max(), Direction::None)))
    {
        quarantine_boundary(quarantine_areas);
    }

    QuarantineEscape() = delete;

    /**
     * Computes whether infection in certain year escaped from quarantine areas
     * and if not, computes and saves minimum distance and direction to quarantine areas
     * for the specified step (simulation year). Aggregates over all quarantine areas.
     */
    void infection_escape_quarantine(
        const IntegerRaster& infected,
        const IntegerRaster& quarantine_areas,
        unsigned simulation_year)
    {
        DistDir min_dist_dir =
            std::make_tuple(std::numeric_limits<double>::max(), Direction::None);
        for (int i = 0; i < height_; i++) {
            for (int j = 0; j < width_; j++) {
                if (!infected(i, j))
                    continue;
                int area = quarantine_areas(i, j);
                if (area == 0) {
                    escape_dist_dirs.at(simulation_year) = std::make_tuple(
                        true, std::make_tuple(std::nan(""), Direction::None));
                    return;
                }
                double dist;
                Direction dir;
                int bindex = boundary_id_idx_map[area];
                std::tie(dist, dir) = closest_direction(i, j, boundaries.at(bindex));
                if (dist < std::get<0>(min_dist_dir)) {
                    min_dist_dir = std::make_tuple(dist, dir);
                }
            }
        }
        escape_dist_dirs.at(simulation_year) = std::make_tuple(false, min_dist_dir);
    }
    /**
     * Computes escape info (if escaped, distance and direction if not escaped)
     * for certain step (simulation year)
     */
    const EscapeDistDir yearly_escape_info(unsigned simulation_year) const
    {
        return escape_dist_dirs.at(simulation_year);
    }
};

/**
 * Reports probability of escaping quarantine based on multiple runs for certain step
 * (simulation year). 1 means 100% probability of escaping.
 */
template<typename IntegerRaster>
double quarantine_escape_probability(
    const std::vector<QuarantineEscape<IntegerRaster>> escape_infos,
    unsigned simulation_year)
{
    bool escape;
    DistDir distdir;
    int escapes = 0;
    for (const auto& item : escape_infos) {
        std::tie(escape, distdir) = item.yearly_escape_info(simulation_year);
        escapes += escape;
    }
    return (double)escapes / escape_infos.size();
}

/**
 * Reports minimum distances to quarantine boundary (bbox)
 * for each run for certain step (simulation year).
 * If in certain runs infection escaped, it reports nan for distance.
 */
template<typename IntegerRaster>
std::vector<double> distance_to_quarantine(
    const std::vector<QuarantineEscape<IntegerRaster>> escape_infos,
    unsigned simulation_year)
{
    bool escape;
    DistDir distdir;
    double dist;
    Direction dir;
    std::vector<double> distances;
    for (const auto& item : escape_infos) {
        std::tie(escape, distdir) = item.yearly_escape_info(simulation_year);
        std::tie(dist, dir) = distdir;
        distances.push_back(dist);
    }

    return distances;
}

/**
 * Writes quarantine escape summary for all years into a string.
 * Uses CSV fomat with commas (year, probability of escape, distances from runs).
 * If escaped in particular run, there is empty value for that distance
 * (...,10,,20,...). It assumes yearly records.
 */
template<typename IntegerRaster>
std::string write_quarantine_escape(
    const std::vector<QuarantineEscape<IntegerRaster>> escape_infos,
    unsigned num_years,
    int start_year)
{
    std::stringstream ss;
    ss << std::setprecision(1) << std::fixed;
    ss << "year,escape_probability";
    for (unsigned i = 0; i < escape_infos.size(); i++)
        ss << ",dist" << i;
    ss << "\n";
    for (unsigned step = 0; step < num_years; step++) {
        ss << start_year + step;
        ss << "," << quarantine_escape_probability(escape_infos, step);
        std::vector<double> dists = distance_to_quarantine(escape_infos, step);
        for (unsigned i = 0; i < dists.size(); i++) {
            if (std::isnan(dists.at(i)))
                ss << ",";
            else
                ss << "," << dists.at(i);
        }
        ss << "\n";
    }
    return ss.str();
}
}  // namespace pops
#endif  // POPS_QUARANTINE_HPP
