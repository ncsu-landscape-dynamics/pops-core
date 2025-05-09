/*
 * PoPS model - network dispersal kernel
 *
 * Copyright (C) 2025 by the authors.
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

#ifndef MULTI_NETWORK_HPP
#define MULTI_NETWORK_HPP

#include "network.hpp"
#include "utils.hpp"

namespace pops {

/**
 * @brief The MultiNetwork class handles one or more networks.
 *
 * @see Network
 */
template<typename RasterIndex>
class MultiNetwork
{
public:
    /**
     * @brief Create MultiNetwork with multiple (empty) networks
     *
     * If a network uses teleport, value of min and max distance is ignored
     * by the individual underlying network.
     *
     * @param bbox Bounding box of the raster grid (in real-world coordinates)
     * @param ew_res East-west resolution of the raster grid
     * @param ns_res North-south resolution of the raster grid
     * @param movements Travel mode for each network (walk, jump, or teleport)
     * @param min_distances Minimum travel distance for each network
     * @param max_distances Maximum travel distance for each network
     */
    MultiNetwork(
        BBox<double> bbox,
        double ew_res,
        double ns_res,
        const std::vector<std::string>& movements,
        const std::vector<double>& min_distances,
        const std::vector<double>& max_distances)
    {
        if (movements.size() != min_distances.size()
            || min_distances.size() != max_distances.size())
            throw std::invalid_argument(std::string(
                "Size of movements (" + std::to_string(movements.size())
                + "), min_distances (" + std::to_string(min_distances.size())
                + "), and max_distances (" + std::to_string(max_distances.size())
                + ") should be the same."));
        for (size_t i = 0; i < movements.size(); ++i) {
            Network<RasterIndex> network(
                bbox, ew_res, ns_res, movements[i], min_distances[i], max_distances[i]);
            networks_.push_back(network);
        }
    }
    /**
     * @brief Load one network from an input stream.
     *
     * The network is identified using an *index*. As for the network order,
     * the indices are the same as indices of the lists provided in the
     * constructor.
     *
     * @param index zero-based index of the network to load the data to
     * @param stream Input stream containing text records for network
     * @param allow_empty True if the loaded network can be empty
     *
     * @see Network::load()
     */
    template<typename InputStream>
    void load(size_t index, InputStream& stream, bool allow_empty = false)
    {
        if (index >= networks_.size()) {
            throw std::out_of_range(
                "Loading to network which was not created (index is "
                + std::to_string(index) + ", but number of networks is "
                + std::to_string(networks_.size()) + ")");
        }
        networks_[index].load(stream, allow_empty);
    }
    /**
     * @brief Move from cell to cell through the network
     *
     * For each network, uses the movement type set in the construtor.
     * For walking and jumping, also the min and max distances will be used.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param generator Random number generator
     * @return Final row and column pair
     *
     * @see Network::move()
     */
    template<typename Generator>
    std::tuple<int, int> move(int row, int col, Generator& generator) const
    {
        auto selected_network =
            pick_network(row, col, generator); /* needs to be eligible */
        return selected_network->move(row, col, generator);
    }
    /**
     * @brief Pick a random network from the contained networks
     *
     * Network is considered selectable only if it has a node at the given cell.
     * Test beforehand using has_node_at() if a node exists.
     *
     * @param row Row index of the cell
     * @param col Column index of the cell
     * @param generator Random number generator
     *
     * @return
     *
     * @see has_node_at()
     */
    template<typename Generator>
    const Network<RasterIndex>*
    pick_network(int row, int col, Generator& generator) const
    {
        std::vector<const Network<RasterIndex>*> selectable_networks;
        selectable_networks.reserve(networks_.size());
        for (const auto& network : networks_) {
            if (network.has_node_at(row, col)) {
                selectable_networks.push_back(&network);
            }
        }
        return pick_random_item(selectable_networks, generator);
    }
    /**
     * @brief Test if the network has a node at a given row and column
     * @param row Row
     * @param col Column
     * @return True if there is at least one node for one network
     *
     * @see Network::has_node_at()
     */
    bool has_node_at(RasterIndex row, RasterIndex col) const
    {
        for (const auto& network : networks_) {
            if (network.has_node_at(row, col)) {
                return true;
            }
        }
        return false;
    }

protected:
    /** Reference to the network */
    std::vector<Network<RasterIndex>> networks_;
    /** Travel distance (cost) distribution */
    std::uniform_real_distribution<double> distance_distribution_;
    /** Step through network instead of traveling between nodes */
    bool teleport_{false};
    /** Snap to nodes when traveling between nodes */
    bool jump_{false};
};

}  // namespace pops

#endif  // MULTI_NETWORK_HPP
