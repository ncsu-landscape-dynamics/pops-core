#ifdef POPS_TEST

/*
 * Simple compilation test for the PoPS spread rate class.
 *
 * Copyright (C) 2018-2019 by the authors.
 *
 * Authors: Anna Petrasova <akratoc gmail com>
 *          Vaclav Petras <wenzeslaus gmail com>
 *
 * This file is part of PoPS.

 * PoPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * PoPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with PoPS. If not, see <https://www.gnu.org/licenses/>.
 */

#include <tuple>
#include <vector>
#include <cmath>

#include <pops/raster.hpp>
#include <pops/spread_rate.hpp>

using namespace pops;

/**
 * Class pretending to be host pool only with methods
 * needed by SpreadRateAction and a method for setting infected raster.
 */
class SpreadRateHostPool
{
public:
    SpreadRateHostPool(
        Raster<int> infected, std::vector<std::vector<int>>& suitable_cells)
        : infected_(infected), suitable_cells_(suitable_cells)
    {}
    void set_infected(Raster<int> infected)
    {
        infected_ = infected;
    }
    int infected_at(int row, int col) const
    {
        return infected_(row, col);
    }
    const std::vector<std::vector<int>>& suitable_cells() const
    {
        return suitable_cells_;
    }

private:
    Raster<int> infected_;
    std::vector<std::vector<int>>& suitable_cells_;
};

int test_spread_rate()
{
    int err = 0;
    Raster<int> infected = {
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0}};

    Raster<int> infected1 = {
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 1, 2, 7, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0}};

    Raster<int> infected2 = {
        {0, 0, 0, 0, 0},
        {0, 0, 2, 0, 0},
        {1, 1, 2, 7, 0},
        {0, 0, 0, 9, 0},
        {0, 0, 0, 0, 0}};

    Raster<int> infected3 = {
        {0, 0, 0, 0, 0},
        {0, 0, 2, 0, 1},
        {1, 1, 0, 7, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0}};

    std::vector<std::vector<int>> suitable_cells = {
        {0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 0}, {1, 1}, {1, 2}, {1, 3},
        {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {3, 0}, {3, 1}, {3, 2},
        {3, 3}, {3, 4}, {4, 0}, {4, 1}, {4, 2}, {4, 3}, {4, 4}};

    SpreadRateHostPool host_pool(infected, suitable_cells);
    SpreadRateAction<SpreadRateHostPool, int> spread_rate(
        host_pool, infected.rows(), infected.cols(), 10, 10, 3);
    host_pool.set_infected(infected1);
    spread_rate.action(host_pool, 0);
    host_pool.set_infected(infected2);
    spread_rate.action(host_pool, 1);
    host_pool.set_infected(infected3);
    spread_rate.action(host_pool, 2);

    double n, s, e, w;
    std::tie(n, s, e, w) = spread_rate.step_rate(0);
    if (!(n == 0 && s == 0 && e == 10 && w == 10)) {
        std::cout << "spread rate for year 1 fails" << std::endl;
        err++;
    }
    std::tie(n, s, e, w) = spread_rate.step_rate(1);
    if (!(n == 10 && s == 10 && e == 0 && w == 10)) {
        std::cout << "spread rate for year 1 fails" << std::endl;
        err++;
    }
    std::tie(n, s, e, w) = spread_rate.step_rate(2);
    if (!(n == 0 && s == -10 && e == 10 && std::isnan(w))) {
        std::cout << "spread rate for year 1 fails" << std::endl;
        err++;
    }

    return err;
}

int main()
{
    int num_errors = 0;

    num_errors += test_spread_rate();
    std::cout << "Spread rate number of errors: " << num_errors << std::endl;
    return num_errors;
}

#endif  // POPS_TEST
