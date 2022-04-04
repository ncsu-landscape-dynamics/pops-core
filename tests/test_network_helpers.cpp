#include <fstream>
#include <regex>
#include <random>

#include <pops/network.hpp>

using namespace pops;

int test_edge_geometry()
{
    int ret = 0;
    using Cell = std::pair<int, int>;
    EdgeGeometry<Cell> geometry;
    geometry.push_back({1, 2});
    geometry.push_back({2, 2});
    geometry.push_back({3, 2});
    geometry.push_back({4, 2});
    geometry.set_cost_per_cell(5);

    if (geometry.cost_per_cell() != 5) {
        std::cerr << "Total cost is wrong: " << geometry.cost_per_cell() << "\n";
        ++ret;
    }
    if (geometry.cost() != 5 * 3) {
        std::cerr << "Total cost is wrong: " << geometry.cost() << "\n";
        ++ret;
    }

    std::vector<std::pair<double, EdgeGeometry<Cell>::size_type>> costs_and_indices = {
        {0.0, 0},
        {1.0, 0},
        {2.4, 0},
        {2.5, 1},
        {4.0, 1},
        {5.0, 1},
        {5.4, 1},
        {5.5, 1},
        {6.0, 1},
        {7.5, 2},
        {10.0, 2},
        {12.0, 2},
        {12.5, 3},
        {15.0, 3},
        {17.0, 3},
        {17.4, 3}};
    for (auto item : costs_and_indices) {
        auto cost = item.first;
        auto expected_index = item.second;
        auto index = geometry.index_from_cost(cost);
        if (index != expected_index) {
            std::cerr << "Got index " << index << " but expected index "
                      << expected_index << " for cost " << cost << "\n";
            ++ret;
        }
    }

    std::vector<std::pair<double, Cell>> costs_and_cells = {
        {0.0, {1, 2}},
        {1.0, {1, 2}},
        {5.0, {2, 2}},
        {7.5, {3, 2}},
        {10.0, {3, 2}},
        {12.4, {3, 2}},
        {12.5, {4, 2}},
        {15.0, {4, 2}},
        {17.4, {4, 2}}};
    for (auto item : costs_and_cells) {
        auto cost = item.first;
        auto expected_cell = item.second;
        auto cell = geometry.cell_by_cost(cost);
        if (cell != expected_cell) {
            std::cerr << "Got cell (" << cell.first << ", " << cell.second
                      << ") but expected cell (" << expected_cell.first << ", "
                      << expected_cell.second << ") for cost " << cost << "\n";
            ++ret;
        }
    }

    EdgeGeometryView<EdgeGeometry<Cell>> view(
        geometry.begin(), geometry.end(), geometry);
    if (view.front() != std::make_pair(1, 2)) {
        std::cerr << "First view item is: (" << view.front().first << ", "
                  << view.front().second << ")\n";
        ++ret;
    }
    if (view.back() != std::make_pair(4, 2)) {
        std::cerr << "Last view item is: (" << view.back().first << ", "
                  << view.back().second << ")\n";
        ++ret;
    }
    if (view.cell_by_cost(0) != std::make_pair(1, 2)) {
        std::cerr << "First view item isby cost: (" << view.cell_by_cost(0).first
                  << ", " << view.cell_by_cost(0).second << ")\n";
        ++ret;
    }
    if (view.cell_by_cost(15) != std::make_pair(4, 2)) {
        std::cerr << "Last view item is by cost: (" << view.cell_by_cost(15).first
                  << ", " << view.cell_by_cost(15).second << ")\n";
        ++ret;
    }

    EdgeGeometryView<EdgeGeometry<Cell>> reversed_view(
        geometry.rbegin(), geometry.rend(), geometry);
    if (reversed_view.front() != std::make_pair(4, 2)) {
        std::cerr << "First reversed_view item is: (" << reversed_view.front().first
                  << ", " << reversed_view.front().second << ")\n";
        ++ret;
    }
    if (reversed_view.back() != std::make_pair(1, 2)) {
        std::cerr << "Last reversed_view item is: (" << reversed_view.back().first
                  << ", " << reversed_view.back().second << ")\n";
        ++ret;
    }
    if (reversed_view.cell_by_cost(0) != std::make_pair(4, 2)) {
        std::cerr << "First reversed_view item isby cost: ("
                  << reversed_view.cell_by_cost(0).first << ", "
                  << reversed_view.cell_by_cost(0).second << ")\n";
        ++ret;
    }
    if (reversed_view.cell_by_cost(15) != std::make_pair(1, 2)) {
        std::cerr << "Last reversed_view item is by cost: ("
                  << reversed_view.cell_by_cost(15).first << ", "
                  << reversed_view.cell_by_cost(15).second << ")\n";
        ++ret;
    }

    return ret;
}

int run_tests()
{
    int ret = 0;

    ret += test_edge_geometry();

    if (ret)
        std::cerr << "Number of errors in the network test: " << ret << "\n";
    return ret;
}

int main()
{
    run_tests();
}
