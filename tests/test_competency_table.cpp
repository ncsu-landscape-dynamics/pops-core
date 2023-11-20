#ifdef POPS_TEST

/*
 * Tests for the competency table functionality.
 *
 * Copyright (C) 2023 by the authors.
 *
 * Authors: Vaclav Petras <wenzeslaus gmail com>
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

#include <pops/competency_table.hpp>
#include <pops/environment.hpp>
#include <pops/raster.hpp>

using namespace pops;

class MockupHostPool : public HostPoolInterface<int>
{
public:
    using Environment = ::Environment<
        Raster<int>,
        Raster<double>,
        Raster<double>::IndexType,
        RandomNumberGeneratorProvider<std::default_random_engine>>;

    Raster<int> hosts;

    int total_hosts_at(int row, int col) const override
    {
        return hosts(row, col);
    }
};

template<typename Table, typename Host>
int check_competency(
    const Table& competency_table,
    int row,
    int col,
    double expected,
    const Host& host_pool,
    std::string name)
{
    double actual_competency = competency_table.competency_at(row, col, &host_pool);
    if (actual_competency != expected) {
        std::cerr << name << ": actual and expected competencies differ at " << row
                  << ", " << col << ": " << actual_competency << " != " << expected
                  << "\n";
        return 1;
    }
    return 0;
}

int test_competency_table_small_example()
{
    int ret = 0;
    MockupHostPool::Environment environment;

    CompetencyTable<MockupHostPool, Raster<double>::IndexType> competency_table(
        environment);

    MockupHostPool host_pool_1;
    MockupHostPool host_pool_2;
    MockupHostPool host_pool_3;

    host_pool_1.hosts = {{0, 0, 0}, {5, 5, 5}, {5, 5, 5}};
    host_pool_2.hosts = {{9, 9, 9}, {0, 0, 0}, {9, 9, 9}};
    host_pool_3.hosts = {{0, 0, 3}, {0, 0, 3}, {0, 0, 3}};

    environment.add_host(&host_pool_1);
    environment.add_host(&host_pool_2);
    environment.add_host(&host_pool_3);

    double expected_competency_1 = 0.10;
    double expected_competency_2 = 0.20;
    double expected_competency_3 = 0.30;
    double expected_competency_4 = 0.40;

    competency_table.add_host_competencies({1, 0, 0}, expected_competency_1);
    competency_table.add_host_competencies({0, 1, 0}, expected_competency_2);
    competency_table.add_host_competencies({1, 1, 0}, expected_competency_3);
    competency_table.add_host_competencies({1, 1, 1}, expected_competency_4);

    ret += check_competency(
        competency_table,
        0,
        0,
        0,
        host_pool_1,
        "test_competency_table_small_example host_pool_1");
    ret += check_competency(
        competency_table,
        0,
        0,
        expected_competency_2,
        host_pool_2,
        "test_competency_table_small_example host_pool_2");
    // We could also do host_pool = nullprt, then it would be competency of the cell.
    ret += check_competency(
        competency_table,
        1,
        1,
        expected_competency_1,
        host_pool_1,
        "test_competency_table_small_example host_pool_1");
    ret += check_competency(
        competency_table,
        1,
        1,
        0,
        host_pool_3,
        "test_competency_table_small_example host_pool_3");
    ret += check_competency(
        competency_table,
        2,
        1,
        expected_competency_3,
        host_pool_1,
        "test_competency_table_small_example host_pool_1");
    ret += check_competency(
        competency_table,
        2,
        2,
        expected_competency_4,
        host_pool_1,
        "test_competency_table_small_example host_pool_1");
    return ret;
}

int main()
{
    int ret = 0;

    ret += test_competency_table_small_example();
    std::cout << "Test of competency table: number of errors: " << ret << std::endl;

    return ret;
}

#endif  // POPS_TEST
