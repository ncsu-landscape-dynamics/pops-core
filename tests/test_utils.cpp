/*
 * Tests of various general utilities
 *
 * Copyright (C) 2025 by the authors.
 *
 * Authors: Vaclav Petras <wenzeslaus gmail com>
 *
 * This file is part of PoPS.
 *
 * PoPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * PoPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PoPS. If not, see <https://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <regex>
#include <random>

#include <pops/model.hpp>

using namespace pops;

#include <random>

/**
 * Test uniform random selection of items
 */
int test_pick_random_item()
{
    int ret{0};
    std::vector<int> numbers{11, 12, 13, 14};
    std::default_random_engine generator(1);
    std::map<int, int> frequencies;
    // While std::map auto-inits to 0,
    // we need to create the records with zero in case they are not created later
    // (because then they would simply not exist in the map and slip the test).
    for (auto number : numbers) {
        frequencies[number] = 0;
    }
    for (int i = 0; i < 1e3; ++i) {
        auto number = pick_random_item(numbers, generator);
        ++frequencies[number];
    }
    std::cout << "test_pick_random_item:\n";
    for (auto const& item : frequencies) {
        std::cout << item.first << ": " << item.second << "x\n";
        int min_expected = 220;
        if (item.second < min_expected) {
            std::cerr
                << "test_pick_random_item: We expect the number to be picked more than "
                << min_expected << " times\n";
            ++ret;
        }
    }
    return ret;
}

/**
 * Test weighted random selection of items
 */
int test_pick_weighted_random_item()
{
    int ret{0};
    std::vector<int> numbers{11, 12, 13, 14};
    // The weights are trivial, but we are using library functions underneath, so here
    // we test only that we are passing values correctly.
    std::vector<double> weights{1, 0, 1, 0};
    std::default_random_engine generator(1);
    std::map<int, int> frequencies;
    // While std::map auto-inits to 0,
    // we need to create the records with zero in case they are not created later
    // (because then they would simply not exist in the map and slip the test).
    for (auto number : numbers) {
        frequencies[number] = 0;
    }
    for (int i = 0; i < 1e3; ++i) {
        auto number = pick_weighted_random_item(numbers, weights, generator);
        ++frequencies[number];
    }
    std::cout << "test_pick_weighted_random_item:\n";
    for (auto const& item : frequencies) {
        std::cout << item.first << ": " << item.second << "x\n";
    }
    int min_expected = 490;
    for (auto number : {11, 13}) {
        if (frequencies[number] < min_expected) {
            std::cerr << "test_pick_weighted_random_item: We expect the number to be "
                         "picked more than "
                      << min_expected << " times\n";
            ++ret;
        }
    }
    for (auto number : {12, 14}) {
        if (frequencies[number] != 0) {
            std::cerr << "test_pick_weighted_random_item: We don't expect the number "
                         "to be picked "
                      << min_expected << " times\n";
            ++ret;
        }
    }
    return ret;
}

int run_tests()
{
    int ret = 0;

    ret += test_pick_random_item();
    ret += test_pick_weighted_random_item();

    if (ret)
        std::cerr << "Number of errors in the utils test: " << ret << "\n";
    return ret;
}

int main()
{
    return run_tests();
}
