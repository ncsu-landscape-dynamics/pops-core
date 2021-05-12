/*
 * PoPS model - random uniform dispersal kernel
 *
 * Copyright (C) 2015-2020 by the authors.
 *
 * Authors: Vaclav Petras (wenzeslaus gmail com)
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_NETWORK_KERNEL_HPP
#define POPS_NETWORK_KERNEL_HPP

#include "kernel_types.hpp"
#include "raster.hpp"
#include "utils.hpp"

#include <set>
#include <random>
#include <string>
#include <sstream>
#include <map>
#include <cmath>
#include <vector>

namespace pops {

template<typename RasterIndex>
class Network
{
public:
    using NodeId = int;
    using Statistics = std::map<std::string, int>;

    Network(BBox<double> bbox, double ew_res, double ns_res, double speed)
        : bbox_(bbox),
          ew_res_(ew_res),
          ns_res_(ns_res),
          cell_travel_time_(((ew_res + ns_res) / 2) / speed)
    {}

    std::pair<RasterIndex, RasterIndex> xy_to_row_col(double x, double y) const
    {
        double col = (x - bbox_.west) / ew_res_;
        double row = (bbox_.north - y) / ns_res_;
        return {std::floor(row), std::floor(col)};
    }

    std::pair<RasterIndex, RasterIndex>
    xy_to_row_col(const std::string& x, const std::string& y) const
    {
        return xy_to_row_col(std::stod(x), std::stod(y));
    }

    NodeId node_id_from_text(const std::string& text) const
    {
        return std::stoi(text);
    }

    bool out_of_bbox(double x, double y) const
    {
        return x > bbox_.east || x < bbox_.west || y > bbox_.north || y < bbox_.south;
    }

    /**
     *
     * std::string node_file, std::string segment_file
     * std::ifstream node_stream{node_file};
     */
    template<typename InputStream>
    void load(
        InputStream& node_stream, InputStream& segment_stream, bool allow_empty = false)
    {
        std::set<NodeId> node_ids;
        load_nodes(node_stream, node_ids);
        if (node_ids.empty()) {
            if (allow_empty)
                return;
            else
                throw std::runtime_error("Network: No nodes within the extend");
        }
        load_segments(segment_stream, node_ids);
    }

    std::vector<NodeId> candidate_nodes_from_node_matrix(NodeId node) const
    {
        std::vector<NodeId> nodes;
        for (const auto& item : node_matrix_) {
            if (item.first == node)
                nodes.push_back(item.second);
            else if (item.second == node)
                nodes.push_back(item.first);
        }
        return nodes;
    }

    template<typename Generator>
    NodeId next_node(NodeId start, Generator& generator) const
    {
        auto nodes = candidate_nodes_from_node_matrix(start);
        return pick_random_node(nodes, generator);
    }

    std::set<NodeId> get_nodes_at(RasterIndex row, RasterIndex col) const
    {
        // return nodes_by_row_col_.at(std::make_pair(row, col));
        auto it = nodes_by_row_col_.find(std::make_pair(row, col));
        if (it != nodes_by_row_col_.end())
            return it->second;
        return std::set<NodeId>();
    }

    template<typename Generator>
    NodeId
    get_random_node_at(RasterIndex row, RasterIndex col, Generator& generator) const
    {
        auto nodes = get_nodes_at(row, col);
        auto num_nodes = nodes.size();
        if (num_nodes == 1) {
            return *nodes.begin();
        }
        else if (num_nodes > 1) {
            return pick_random_node(nodes, generator);
        }
        else {
            throw std::invalid_argument("No nodes at a given row and column");
        }
    }

    bool has_node_at(RasterIndex row, RasterIndex col) const
    {
        // Replace by contains in C++20.
        // return nodes_by_row_col_.count(std::make_pair(row, col)) > 0;
        auto it = nodes_by_row_col_.find(std::make_pair(row, col));
        if (it == nodes_by_row_col_.end())
            return false;
        return true;
    }

    template<typename Generator>
    std::tuple<int, int> travel(
        RasterIndex start_row,
        RasterIndex start_col,
        double time,
        Generator& generator) const
    {
        // We assume there is a node here, i.e., that we are made decision
        // to use this kernel knowing there is a node.
        auto node_id = get_random_node_at(start_row, start_col, generator);
        while (time >= 0) {
            auto next_node_id = next_node(node_id, generator);
            auto segment = segments_by_nodes_.at(std::make_pair(node_id, next_node_id));
            // nodes may need special handling
            for (auto cell : segment) {
                time -= cell_travel_time_;
                if (time <= 0) {
                    return cell;
                    // Given the while condition, this subsequently ends the while loop
                    // as well.
                    // break;
                }
            }
            node_id = next_node_id;
        }
        throw std::invalid_argument("Time must be greater than or equal to zero");
    }

    Statistics collect_statistics() const
    {
        std::map<std::string, int> stats;
        std::set<NodeId> node_ids;
        for (auto item : nodes_by_row_col_) {
            for (auto node_id : item.second) {
                node_ids.insert(node_id);
            }
        }
        stats["num_nodes"] = node_ids.size();
        stats["num_segments"] = node_matrix_.size();
        std::set<NodeId> nodes_with_segments;
        for (const auto& item : node_matrix_) {
            nodes_with_segments.insert(item.first);
            nodes_with_segments.insert(item.second);
        }
        stats["num_nodes_with_segments"] = nodes_with_segments.size();
        int num_standalone_nodes = 0;
        for (NodeId node_id : node_ids) {
            if (nodes_with_segments.find(node_id) == nodes_with_segments.end()) {
                stats["standalone_node_" + std::to_string(++num_standalone_nodes)] =
                    node_id;
            }
        }
        stats["num_standalone_nodes"] = num_standalone_nodes;
        RasterIndex min_row;
        RasterIndex min_col;
        RasterIndex max_row;
        RasterIndex max_col;
        std::tie(min_row, min_col) = xy_to_row_col(bbox_.west, bbox_.north);
        std::tie(max_row, max_col) = xy_to_row_col(bbox_.east, bbox_.south);
        stats["min_row"] = min_row;
        stats["min_col"] = min_col;
        stats["max_row"] = max_row;
        stats["max_col"] = max_col;
        return stats;
    }

    template<typename OutputStream>
    void dump_yaml(OutputStream& stream) const
    {
        stream << "network:\n";
        stream << "  statistics:\n";
        auto stats = collect_statistics();
        for (const auto& item : stats) {
            stream << "    " << item.first << ": " << item.second << "\n";
        }
        stream << "  edges:\n";
        for (const auto& item : node_matrix_) {
            stream << "    - [" << item.first << ", " << item.second << "]\n";
        }
        stream << "  nodes:\n";
        for (const auto& item : nodes_by_row_col_) {
            auto row = item.first.first;
            auto col = item.first.second;
            for (const auto& node_id : item.second) {
                stream << "    - id: " << node_id << "\n";
                stream << "      row: " << row << "\n";
                stream << "      col: " << col << "\n";
            }
        }
    }

protected:
    template<typename InputStream>
    void load_nodes(InputStream& stream, std::set<NodeId>& node_ids)
    {
        std::string line;
        while (std::getline(stream, line)) {
            std::istringstream line_stream{line};
            char delimeter{','};
            std::string node_text;
            std::getline(line_stream, node_text, delimeter);
            std::string x_coord_text;
            std::getline(line_stream, x_coord_text, delimeter);
            std::string y_coord_text;
            std::getline(line_stream, y_coord_text, delimeter);
            RasterIndex row;
            RasterIndex col;
            double x = std::stod(x_coord_text);
            double y = std::stod(y_coord_text);
            // Cut to extend
            if (out_of_bbox(x, y))
                continue;
            std::tie(row, col) = xy_to_row_col(x_coord_text, y_coord_text);
            NodeId node_id = node_id_from_text(node_text);
            if (node_id < 1) {
                std::runtime_error("Node ID must be greater than zero");
            }
            nodes_by_row_col_[std::make_pair(row, col)].insert(node_id);
            node_ids.insert(node_id);
        }
    }

    template<typename InputStream>
    void load_segments(InputStream& stream, const std::set<NodeId>& node_ids)
    {
        std::string line;
        while (std::getline(stream, line)) {
            std::istringstream line_stream{line};
            char delimeter{','};
            std::string node_1_text;
            std::getline(line_stream, node_1_text, delimeter);
            std::string node_2_text;
            std::getline(line_stream, node_2_text, delimeter);
            auto node_1_id = node_id_from_text(node_1_text);
            auto node_2_id = node_id_from_text(node_2_text);
            // If either end nodes of the segment is not in the extent, skip it.
            // TODO: Part of the segment may still be out, so that needs to be checked.
            // Replace by contains for C++20.
            if (node_ids.find(node_1_id) == node_ids.end()
                || node_ids.find(node_2_id) == node_ids.end())
                continue;
            if (node_1_id == node_2_id) {
                std::runtime_error(
                    std::string("Segment cannot begin and end with the same node: ")
                    + node_1_text + " " + node_2_text);
            }
            std::string segment_text;
            std::getline(line_stream, segment_text, delimeter);
            // We don't know which way the nodes are ordered, so instead of checking
            // the order, we create a symmetric matrix since we allocated the memory
            // anyway.
            node_matrix_.emplace(node_1_id, node_2_id);
            std::istringstream segment_stream{segment_text};
            char in_cell_delimeter{';'};
            std::string x_coord_text;
            std::string y_coord_text;
            Segment segment;
            while (std::getline(segment_stream, x_coord_text, in_cell_delimeter)
                   && std::getline(segment_stream, y_coord_text, in_cell_delimeter)) {
                // The same cell is possibly repeated if raster resolution is lower than
                // the detail ("resolution") of each segment, so we check if the last
                // added coordinate pair converted to same cell.
                auto new_point = xy_to_row_col(x_coord_text, y_coord_text);
                if (segment.empty() || segment.back() != new_point)
                    segment.emplace_back(new_point);
            }
            segments_by_nodes_.emplace(
                std::make_pair(node_1_id, node_2_id), std::move(segment));
        }
    }

    template<typename Container, typename Generator>
    static NodeId pick_random_node(Container nodes, Generator& generator)
    {
        auto num_nodes = nodes.size();  // Replace by std::size in C++17.
        std::uniform_int_distribution<size_t> dist(0, num_nodes - 1);
        auto index = dist(generator);
        // Our lists are expected to be short, so this is expected to be fast for
        // both sets and vectors.
        return *std::next(nodes.begin(), index);
    }

    using NodeMatrix = std::set<std::pair<NodeId, NodeId>>;
    using Segment = std::vector<std::pair<RasterIndex, RasterIndex>>;
    using SegmentsByNodes = std::map<std::pair<NodeId, NodeId>, Segment>;

    BBox<double> bbox_;
    double ew_res_;
    double ns_res_;
    double cell_travel_time_;
    std::map<std::pair<RasterIndex, RasterIndex>, std::set<NodeId>> nodes_by_row_col_;
    NodeMatrix node_matrix_;
    SegmentsByNodes segments_by_nodes_;
};

/*! Dispersal kernel for random uniform dispersal over the whole
 * landscape
 *
 * This class is a good example of how to write a kernel and
 * it is useful for testing due to its simplicity. It tends to generate
 * a lot of spread because it quickly spreads over the landscape.
 * However, it may work as a good starting point for cases where no
 * theory about the spread is available.
 */
class NetworkDispersalKernel
{
public:
    // Other kernels don't have it as a template paramater now, so we just define it
    // to have it for the network definition.
    using RasterIndex = int;

    NetworkDispersalKernel(
        Network<RasterIndex>& network, double min_time, double max_time)
        : network_(network), time_distribution_(min_time, max_time)
    {}

    /*! \copybrief RadialDispersalKernel::operator()()
     *
     */
    template<typename Generator>
    std::tuple<int, int> operator()(Generator& generator, int row, int col)
    {
        double time = time_distribution_(generator);
        std::tie(row, col) = network_.travel(row, col, time, generator);

        return std::make_tuple(row, col);
    }

    bool is_cell_eligible(int row, int col)
    {
        return network_.has_node_at(row, col);
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        return type == DispersalKernelType::Network;
    }

protected:
    Network<RasterIndex>& network_;
    std::uniform_real_distribution<double> time_distribution_;
};

}  // namespace pops

#endif  // POPS_NETWORK_KERNEL_HPP
