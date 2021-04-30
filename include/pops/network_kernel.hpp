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
#include <vector>

namespace pops {

template<typename RasterIndex>
class Network
{
public:
    using NodeId = int;

    Network(BBox<double> bbox, double ew_res, double ns_res)
        : bbox_(bbox), ew_res_(ew_res), ns_res_(ns_res)
    {}

    std::pair<RasterIndex, RasterIndex> xy_to_row_col(double x, double y)
    {
        // TODO: implement
        return {100 * y, 100 * x};
    }

    std::pair<RasterIndex, RasterIndex>
    xy_to_row_col(const std::string& x, const std::string& y)
    {
        return xy_to_row_col(std::stod(x), std::stod(y));
    }

    NodeId node_id_from_text(const std::string& text)
    {
        return std::stoi(text);
    }

    bool out_of_bbox(double x, double y)
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
        std::string line;
        std::set<NodeId> node_ids;
        std::map<NodeId, std::pair<double, double>> node_coords;
        while (std::getline(node_stream, line)) {
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
        if (node_ids.empty()) {
            if (allow_empty)
                return;
            else
                throw std::runtime_error("Network: No nodes within the extend");
        }
        while (std::getline(segment_stream, line)) {
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
            while (std::getline(segment_stream, x_coord_text, in_cell_delimeter)
                   && std::getline(segment_stream, y_coord_text, in_cell_delimeter)) {
            }
        }
    }

    std::vector<NodeId> candidate_nodes_from_node_matrix(NodeId node)
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
    NodeId next_node(NodeId start, Generator& generator)
    {
        auto nodes = candidate_nodes_from_node_matrix(start);
        auto num_nodes = nodes.size();
        std::uniform_int_distribution<size_t> dist(0, num_nodes - 1);
        auto index = dist(generator);
        return nodes[index];
    }

    std::set<NodeId> get_nodes_at(RasterIndex row, RasterIndex col) const
    {
        // return nodes_by_row_col_.at(std::make_pair(row, col));
        auto it = nodes_by_row_col_.find(std::make_pair(row, col));
        if (it != nodes_by_row_col_.end())
            return it->second;
        return std::set<NodeId>();
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
    std::tuple<int, int> time_to_row_col(
        RasterIndex start_row, RasterIndex start_col, double time, Generator& generator)
    {
        // The trouble is when there is no node here, but we are already made decision
        // to use this kernel.
        auto node_id = get_node(start_row, start_col);
        while (time > 0) {
            auto next_node_id = next_node(node_id, generator);
            auto segment = get_segment(node_id, next_node_id);
            // nodes may need special handling
            for (auto cell : segment) {
                time -= get_travel_time(cell);
                if (time <= 0) {
                    return get_row_col(cell);
                    // Given the while condition, this subsequently ends the while loop
                    // as well.
                    // break;
                }
            }
            node_id = next_node_id;
        }
    }

    std::map<std::string, int> collect_stats()
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
        return stats;
    }

protected:
    using NodeMatrix = std::set<std::pair<NodeId, NodeId>>;

    BBox<double> bbox_;
    double ew_res_;
    double ns_res_;
    std::map<std::pair<RasterIndex, RasterIndex>, std::set<NodeId>> nodes_by_row_col_;
    NodeMatrix node_matrix_;
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
        std::tie(row, col) = network_.time_to_row_col(row, col, time, generator);

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
