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

template<typename RasterIndex>
class MultiNetwork {
public:
    MultiNetwork(
        BBox<double> bbox, double ew_res, double ns_res,
        const std::vector<std::string>& movements,
        const std::vector<double>& min_distances,
        const std::vector<double>& max_distances)
    {
        if (movements.size() != min_distances.size() || min_distances.size() != max_distances.size())
            throw std::invalid_argument(std::string("Size of movements (" + std::to_string(movements.size()) + "), min_distances (" + std::to_string(min_distances.size()) + "), and max_distances (" + std::to_string(max_distances.size()) + ") should be the same."));
        for (size_t i = 0; i<movements.size(); ++i) {
            Network<RasterIndex> network(bbox, ew_res, ns_res, movements[i], min_distances[i], max_distances[i]);
            networks_.push_back(network);
        }
    }
    /**
     * @brief load
     * @param stream
     * @param allow_empty
     *
     * @see Network::load()
     */
    template<typename InputStream>
    void load(size_t index, InputStream& stream, bool allow_empty = false)
    {
        if (index >= networks_.size()) {
            throw std::out_of_range("Loading to network which was not created (index is " + std::to_string(index) + ", but number of networks is " + std::to_string(networks_.size()) + ")");
        }
        networks_[index].load(stream, allow_empty);
    }
    template<typename Generator>
    std::tuple<int, int> move(int row, int col, Generator& generator) const
    {
        auto selected_network = pick_network(row, col, generator);  /* needs to be eligible */
        return selected_network->move(row, col, generator);
    }
    template<typename Generator>
    const Network<RasterIndex>* pick_network(int row, int col, Generator& generator) const {
        std::vector<const Network<RasterIndex>*> selectable_networks;
        selectable_networks.reserve(networks_.size());
        for (const auto& network: networks_) {
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
        for (const auto& network: networks_) {
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
