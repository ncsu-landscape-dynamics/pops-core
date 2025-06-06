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
#include <limits>
#include <type_traits>
#include <sstream>
#include <iomanip>

#include "utils.hpp"

namespace pops {

std::ostream& operator<<(std::ostream& os, const Direction& obj)
{
    os << static_cast<std::underlying_type<Direction>::type>(obj);
    return os;
}

typedef std::tuple<double, Direction> DistDir;
typedef std::tuple<bool, DistDir> EscapeDistDir;
typedef std::vector<EscapeDistDir> EscapeDistDirs;
typedef std::map<Direction, bool> Directions;

/*! Get a map with direction and bool parsed from string with letters N S E W
 *  separated by delimeter. Empty string means all directions are true.
 *
 * Throws an std::invalid_argument exception if the value was not
 * found or is not supported (which is the same thing).
 */
Directions directions_from_string(const std::string& text, char delimiter = ',')
{
    Directions directions;
    for (const auto& direction :
         {Direction::N, Direction::E, Direction::S, Direction::W})
        directions[direction] = text.empty() ? true : false;

    if (text.empty())
        return directions;
    std::string direction;
    std::istringstream directionStream(text);
    while (std::getline(directionStream, direction, delimiter)) {
        std::map<std::string, Direction> mapping{
            {"N", Direction::N},
            {"E", Direction::E},
            {"S", Direction::S},
            {"W", Direction::W}};
        try {
            directions[mapping.at(direction)] = true;
        }
        catch (const std::out_of_range&) {
            throw std::invalid_argument(
                "directions_from_string: Invalid value '" + direction + "' provided");
        }
    }
    return directions;
}

/**
 * Class storing and computing quarantine escape metrics for one simulation.
 */
template<typename IntegerRaster, typename RasterIndex = int>
class QuarantineEscapeAction
{
private:
    RasterIndex width_;
    RasterIndex height_;
    // the west-east resolution of the pixel
    double west_east_resolution_;
    // the north-south resolution of the pixel
    double north_south_resolution_;
    unsigned num_steps;
    std::vector<BBoxInt> boundaries;
    Directions directions_;
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
        if (directions_.at(Direction::N)
            && (i - n) * north_south_resolution_ < mindist) {
            mindist = std::lround((i - n) * north_south_resolution_);
            closest = std::make_tuple(mindist, Direction::N);
        }
        if (directions_.at(Direction::S)
            && (s - i) * north_south_resolution_ < mindist) {
            mindist = std::lround((s - i) * north_south_resolution_);
            closest = std::make_tuple(mindist, Direction::S);
        }
        if (directions_.at(Direction::E) && (e - j) * west_east_resolution_ < mindist) {
            mindist = std::lround((e - j) * west_east_resolution_);
            closest = std::make_tuple(mindist, Direction::E);
        }
        if (directions_.at(Direction::W) && (j - w) * west_east_resolution_ < mindist) {
            mindist = std::lround((j - w) * west_east_resolution_);
            closest = std::make_tuple(mindist, Direction::W);
        }
        return closest;
    }

public:
    QuarantineEscapeAction(
        const IntegerRaster& quarantine_areas,
        double ew_res,
        double ns_res,
        unsigned num_steps,
        std::string directions = "")
        : width_(quarantine_areas.cols()),
          height_(quarantine_areas.rows()),
          west_east_resolution_(ew_res),
          north_south_resolution_(ns_res),
          num_steps(num_steps),
          directions_(directions_from_string(directions)),
          escape_dist_dirs(
              num_steps,
              std::make_tuple(
                  false,
                  std::make_tuple(std::numeric_limits<double>::max(), Direction::None)))
    {
        quarantine_boundary(quarantine_areas);
    }

    QuarantineEscapeAction() = delete;

    /**
     * Computes whether infection in certain step escaped from quarantine areas
     * and if not, computes and saves minimum distance and direction to quarantine areas
     * for the specified step. Aggregates over all quarantine areas.
     */
    template<typename Hosts>
    void
    action(const Hosts& hosts, const IntegerRaster& quarantine_areas, unsigned step)
    {

        DistDir min_dist_dir =
            std::make_tuple(std::numeric_limits<double>::max(), Direction::None);
        for (auto indices : hosts.suitable_cells()) {
            int i = indices[0];
            int j = indices[1];
            if (!hosts.infected_at(i, j))
                continue;
            int area = quarantine_areas(i, j);
            if (area == 0) {
                escape_dist_dirs.at(step) = std::make_tuple(
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
        escape_dist_dirs.at(step) = std::make_tuple(false, min_dist_dir);
    }
    /**
     * Computes escape info (if escaped, distance and direction if not escaped)
     * for certain action step.
     */
    const EscapeDistDir escape_info(unsigned step) const
    {
        return escape_dist_dirs.at(step);
    }

    /**
     * Returns true if infection escaped the quarantine boundary.
     */
    bool escaped(unsigned step) const
    {
        return std::get<0>(escape_dist_dirs.at(step));
    }

    /**
     * Returns minimum distance to quarantine boundary bbox.
     * If infection already escaped, returns NaN.
     * Aggregated over all quarantine areas.
     */
    double distance(unsigned step) const
    {
        auto dist_dir = std::get<1>(escape_dist_dirs.at(step));
        return std::get<0>(dist_dir);
    }

    /**
     * Returns the direction (N, S, E, W, None) of the minimum distance to quarantine
     * boundary bbox. Returns None if infection already escaped.
     */
    Direction direction(unsigned step) const
    {
        auto dist_dir = std::get<1>(escape_dist_dirs.at(step));
        return std::get<1>(dist_dir);
    }
};

/**
 * Reports probability of escaping quarantine based on multiple runs for certain step.
 * 1 means 100% probability of escaping.
 */
template<typename IntegerRaster>
double quarantine_escape_probability(
    const std::vector<QuarantineEscapeAction<IntegerRaster>> escape_infos,
    unsigned step)
{
    bool escape;
    DistDir distdir;
    int escapes = 0;
    for (const auto& item : escape_infos) {
        std::tie(escape, distdir) = item.escape_info(step);
        escapes += escape;
    }
    return (double)escapes / escape_infos.size();
}

/**
 * Reports minimum distances to quarantine boundary (bbox) and associated distances
 * for each run for certain step.
 * If in certain runs infection escaped, it reports nan for distance and None for
 * direction.
 */
template<typename IntegerRaster>
std::vector<DistDir> distance_direction_to_quarantine(
    const std::vector<QuarantineEscapeAction<IntegerRaster>> escape_infos,
    unsigned step)
{
    bool escape;
    DistDir distdir;
    std::vector<DistDir> distances_directions;
    for (const auto& item : escape_infos) {
        std::tie(escape, distdir) = item.escape_info(step);
        distances_directions.push_back(distdir);
    }

    return distances_directions;
}

/**
 * Writes quarantine escape summary for all steps into a string.
 * Uses CSV format with commas (step, probability of escape, distance, direction as
 * azimuth from runs). E.g. "0,0.4,1000,0,2000,90,1500,0" If escaped in particular run,
 * there is empty value for that distance and direction
 * (...,1000,0,,,2000,90,...).
 */
template<typename IntegerRaster>
std::string write_quarantine_escape(
    const std::vector<QuarantineEscapeAction<IntegerRaster>> escape_infos,
    unsigned num_steps)
{
    std::stringstream ss;
    ss << std::setprecision(1) << std::fixed;
    ss << "step,escape_probability";
    for (unsigned i = 0; i < escape_infos.size(); i++)
        ss << ",dist" << i << ",dir" << i;
    ss << "\n";
    for (unsigned step = 0; step < num_steps; step++) {
        ss << step;
        ss << "," << quarantine_escape_probability(escape_infos, step);
        std::vector<DistDir> dists_dirs =
            distance_direction_to_quarantine(escape_infos, step);
        double dist;
        Direction dir;
        for (unsigned i = 0; i < dists_dirs.size(); i++) {
            std::tie(dist, dir) = dists_dirs.at(i);
            if (std::isnan(dist))
                ss << ",,";
            else
                ss << "," << dist << "," << dir;
        }
        ss << "\n";
    }
    return ss.str();
}
}  // namespace pops
#endif  // POPS_QUARANTINE_HPP
