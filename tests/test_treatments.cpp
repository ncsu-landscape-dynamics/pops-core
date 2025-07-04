#ifdef POPS_TEST

/*
 * Simple compilation test for the PoPS Treatments class.
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

#include <vector>
#include <pops/raster.hpp>
#include <pops/treatments.hpp>
#include <pops/date.hpp>
#include <pops/scheduling.hpp>
#include <pops/environment.hpp>

using namespace pops;

using StandardSingleHostPool = HostPool<
    Raster<int>,
    Raster<double>,
    int,
    RandomNumberGeneratorProvider<std::default_random_engine>>;

using TestEnvironment = Environment<
    Raster<int>,
    Raster<double>,
    int,
    RandomNumberGeneratorProvider<std::default_random_engine>>;

int test_application_ratio()
{
    int num_errors = 0;
    Scheduler scheduler(Date(2020, 1, 1), Date(2020, 12, 31), StepUnit::Month, 1);

    TestEnvironment environment;

    Raster<double> tr1 = {{1, 0.5}, {0.75, 0}};
    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    Raster<int> zeros(infected.rows(), infected.cols(), 0);
    auto total_hosts = infected + susceptible + resistant;
    std::vector<Raster<int>> exposed;
    std::vector<Raster<int>> mortality_tracker(1, infected);

    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};

    StandardSingleHostPool host_pool(
        ModelType::SusceptibleInfected,
        bool(mortality_tracker.size()),
        susceptible,
        exposed,
        0,
        infected,
        zeros,
        resistant,
        mortality_tracker,
        zeros,
        total_hosts,
        environment,
        false,
        0,
        false,
        0,
        infected.rows(),
        infected.cols(),
        suitable_cells);

    Treatments<StandardSingleHostPool, Raster<double>> treatments(scheduler);
    treatments.add_treatment(tr1, Date(2020, 1, 1), 0, TreatmentApplication::Ratio);
    treatments.manage(0, host_pool);

    Raster<int> treated = {{0, 3}, {5, 42}};
    Raster<int> inf_treated = {{0, 2}, {4, 40}};
    auto th_treated = treated + inf_treated + resistant;
    if (!(susceptible == treated && infected == inf_treated
          && total_hosts == th_treated)) {
        std::cerr << "Treatment with ratio app does not work\n";
        std::cerr << susceptible << infected << total_hosts;
        num_errors++;
    }
    return num_errors;
}

/**
 * @brief Test treatment without active mortality
 * @return Number of errors encountered
 */
int test_application_ratio_without_mortality()
{
    int num_errors = 0;
    Scheduler scheduler(Date(2020, 1, 1), Date(2020, 12, 31), StepUnit::Month, 1);

    TestEnvironment environment;

    Raster<double> tr1 = {{1, 0.5}, {0.75, 0}};
    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    Raster<int> zeros(infected.rows(), infected.cols(), 0);
    auto total_hosts = infected + susceptible + resistant;
    std::vector<Raster<int>> exposed;
    std::vector<Raster<int>> mortality_tracker;

    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};

    StandardSingleHostPool host_pool(
        ModelType::SusceptibleInfected,
        false,
        susceptible,
        exposed,
        0,
        infected,
        zeros,
        resistant,
        mortality_tracker,
        zeros,
        total_hosts,
        environment,
        false,
        0,
        false,
        0,
        infected.rows(),
        infected.cols(),
        suitable_cells);

    // First, test that host pool works for the case without mortality.
    for (int row = 0; row < infected.rows(); ++row) {
        for (int col = 0; col < infected.cols(); ++col) {
            auto mortality_groups = host_pool.mortality_by_group_at(row, col);
            if (mortality_groups.size() != 1) {
                std::cerr << "Expected a single mortality group from host pool but got "
                          << mortality_groups.size() << "\n";
                num_errors++;
            }
            if (mortality_groups[0] != host_pool.infected_at(row, col)) {
                std::cerr
                    << "Host pool does not work as expected: "
                    << "The single mortality group is diferent from total infected ("
                    << mortality_groups[0] << " != " << host_pool.infected_at(row, col)
                    << ")\n";
                num_errors++;
            }
        }
    }

    Treatments<StandardSingleHostPool, Raster<double>> treatments(scheduler);
    treatments.add_treatment(tr1, Date(2020, 1, 1), 0, TreatmentApplication::Ratio);
    treatments.manage(0, host_pool);

    Raster<int> treated = {{0, 3}, {5, 42}};
    Raster<int> inf_treated = {{0, 2}, {4, 40}};
    auto th_treated = treated + inf_treated + resistant;
    if (!(susceptible == treated && infected == inf_treated
          && total_hosts == th_treated)) {
        std::cerr << "Treatment with ratio app does not work without mortality\n";
        std::cerr << susceptible << infected << total_hosts;
        num_errors++;
    }
    return num_errors;
}

int test_application_all_inf()
{
    int num_errors = 0;
    Scheduler scheduler(Date(2020, 1, 1), Date(2020, 12, 31), StepUnit::Month, 1);

    TestEnvironment environment;

    Raster<double> tr1 = {{1, 0.5}, {0.75, 0}};
    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    Raster<int> zeros(infected.rows(), infected.cols(), 0);
    auto total_hosts = infected + susceptible + resistant;
    std::vector<Raster<int>> exposed;
    std::vector<Raster<int>> mortality_tracker(1, infected);

    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};

    StandardSingleHostPool host_pool(
        ModelType::SusceptibleInfected,
        bool(mortality_tracker.size()),
        susceptible,
        exposed,
        0,
        infected,
        zeros,
        resistant,
        mortality_tracker,
        zeros,
        total_hosts,
        environment,
        false,
        0,
        false,
        0,
        infected.rows(),
        infected.cols(),
        suitable_cells);

    Treatments<StandardSingleHostPool, Raster<double>> treatments(scheduler);
    treatments.add_treatment(
        tr1, Date(2020, 1, 1), 0, TreatmentApplication::AllInfectedInCell);
    treatments.manage(0, host_pool);

    Raster<int> treated = {{0, 3}, {5, 42}};
    Raster<int> inf_treated = {{0, 0}, {0, 40}};
    auto th_treated = treated + inf_treated + resistant;
    if (!(susceptible == treated && infected == inf_treated
          && total_hosts == th_treated)) {
        std::cerr << "Treatment with AllInfectedInCell app does not work\n";
        std::cerr << susceptible << infected << total_hosts;
        num_errors++;
    }
    return num_errors;
}

int test_application_ratio_pesticide()
{
    int num_errors = 0;
    Scheduler scheduler(Date(2020, 1, 1), Date(2020, 12, 31), StepUnit::Day, 7);
    unsigned n = scheduler.schedule_action_date(Date(2020, 1, 1));

    TestEnvironment environment;

    Raster<double> tr1 = {{1, 0.5}, {0.75, 0}};
    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    Raster<int> zeros(infected.rows(), infected.cols(), 0);
    auto total_hosts = infected + susceptible + resistant;
    std::vector<Raster<int>> exposed;
    std::vector<Raster<int>> mortality_tracker(1, infected);

    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};

    StandardSingleHostPool host_pool(
        ModelType::SusceptibleInfected,
        bool(mortality_tracker.size()),
        susceptible,
        exposed,
        0,
        infected,
        zeros,
        resistant,
        mortality_tracker,
        zeros,
        total_hosts,
        environment,
        false,
        0,
        false,
        0,
        infected.rows(),
        infected.cols(),
        suitable_cells);

    Treatments<StandardSingleHostPool, Raster<double>> treatments(scheduler);
    treatments.add_treatment(tr1, Date(2020, 5, 1), 7, TreatmentApplication::Ratio);
    treatments.manage(n, host_pool);

    Raster<int> treated = {{10, 6}, {20, 42}};
    Raster<int> inf_treated = {{1, 4}, {16, 40}};
    Raster<int> resist = {{0, 0}, {0, 0}};
    auto th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr << "Treatment with pesticide and ratio app does not work - block 1\n";
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }
    n = scheduler.schedule_action_date(Date(2020, 5, 3));
    treatments.manage(n, host_pool);

    treated = {{0, 3}, {5, 42}};
    inf_treated = {{0, 2}, {4, 40}};
    resist = {{11, 5}, {27, 0}};
    th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr << "Treatment with pesticide and ratio app does not work - block 2\n";
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }
    n = scheduler.schedule_action_date(Date(2020, 5, 8));
    treatments.manage(n, host_pool);

    treated = {{11, 8}, {32, 42}};
    resist = {{0, 0}, {0, 0}};
    th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr << "Treatment with pesticide and ratio app does not work"
                  << std::endl;
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }
    return num_errors;
}

int test_application_all_inf_pesticide()
{
    int num_errors = 0;
    Scheduler scheduler(Date(2020, 1, 1), Date(2020, 12, 31), StepUnit::Day, 7);

    TestEnvironment environment;

    Raster<double> tr1 = {{1, 0.5}, {0.75, 0}};
    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    Raster<int> zeros(infected.rows(), infected.cols(), 0);
    auto total_hosts = infected + susceptible + resistant;
    std::vector<Raster<int>> exposed;
    std::vector<Raster<int>> mortality_tracker(1, infected);

    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};

    StandardSingleHostPool host_pool(
        ModelType::SusceptibleInfected,
        bool(mortality_tracker.size()),
        susceptible,
        exposed,
        0,
        infected,
        zeros,
        resistant,
        mortality_tracker,
        zeros,
        total_hosts,
        environment,
        false,
        0,
        false,
        0,
        infected.rows(),
        infected.cols(),
        suitable_cells);

    Treatments<StandardSingleHostPool, Raster<double>> treatments(scheduler);
    treatments.add_treatment(
        tr1, Date(2020, 5, 1), 7, TreatmentApplication::AllInfectedInCell);
    unsigned n = scheduler.schedule_action_date(Date(2020, 1, 1));
    treatments.manage(n, host_pool);

    Raster<int> treated = {{10, 6}, {20, 42}};
    Raster<int> inf_treated = {{1, 4}, {16, 40}};
    Raster<int> resist = {{0, 0}, {0, 0}};
    auto th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr << "Treatment with pesticide and AllInfectedInCell app "
                     "does not work - block 1\n";
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }
    n = scheduler.schedule_action_date(Date(2020, 5, 3));
    treatments.manage(n, host_pool);

    treated = {{0, 3}, {5, 42}};
    inf_treated = {{0, 0}, {0, 40}};
    resist = {{11, 7}, {31, 0}};
    th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr << "Treatment with pesticide and AllInfectedInCell app "
                     "does not work - block 2\n";
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }
    n = scheduler.schedule_action_date(Date(2020, 5, 8));
    treatments.manage(n, host_pool);

    treated = {{11, 10}, {36, 42}};
    resist = {{0, 0}, {0, 0}};
    th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr << "Treatment with pesticide and AllInfectedInCell app does not work"
                  << std::endl;
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }
    return num_errors;
}

int test_combination()
{
    int num_errors = 0;
    Scheduler scheduler(Date(2020, 1, 1), Date(2020, 12, 31), StepUnit::Day, 7);
    TestEnvironment environment;

    Raster<double> tr1 = {{1, 0.5}, {0.75, 0}};
    Raster<double> tr2 = {{1, 1}, {1, 1}};

    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    Raster<int> zeros(infected.rows(), infected.cols(), 0);
    auto total_hosts = infected + susceptible + resistant;
    std::vector<Raster<int>> exposed;
    std::vector<Raster<int>> mortality_tracker(1, infected);

    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};

    StandardSingleHostPool host_pool(
        ModelType::SusceptibleInfected,
        bool(mortality_tracker.size()),
        susceptible,
        exposed,
        0,
        infected,
        zeros,
        resistant,
        mortality_tracker,
        zeros,
        total_hosts,
        environment,
        false,
        0,
        false,
        0,
        infected.rows(),
        infected.cols(),
        suitable_cells);

    Treatments<StandardSingleHostPool, Raster<double>> treatments(scheduler);
    treatments.add_treatment(tr1, Date(2020, 5, 1), 0, TreatmentApplication::Ratio);
    treatments.add_treatment(tr2, Date(2020, 6, 1), 7, TreatmentApplication::Ratio);
    unsigned n = scheduler.schedule_action_date(Date(2020, 1, 1));
    treatments.manage(n, host_pool);

    Raster<int> treated = {{10, 6}, {20, 42}};
    Raster<int> inf_treated = {{1, 4}, {16, 40}};
    Raster<int> resist = {{0, 0}, {0, 0}};
    auto th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr << "Combination of treatments does not work - block 1\n";
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }
    n = scheduler.schedule_action_date(Date(2020, 5, 3));
    treatments.manage(n, host_pool);

    treated = {{0, 3}, {5, 42}};
    inf_treated = {{0, 2}, {4, 40}};
    resist = {{0, 0}, {0, 0}};
    th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr << "Combination of treatments does not work - block 2\n";
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }
    n = scheduler.schedule_action_date(Date(2020, 6, 2));
    treatments.manage(n, host_pool);

    treated = {{0, 0}, {0, 0}};
    inf_treated = {{0, 0}, {0, 0}};
    resist = {{0, 5}, {9, 82}};
    th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr << "Combination of treatments does not work - block 3\n";
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }
    n = scheduler.schedule_action_date(Date(2020, 6, 8));
    treatments.manage(n, host_pool);

    treated = {{0, 5}, {9, 82}};
    inf_treated = {{0, 0}, {0, 0}};
    resist = {{0, 0}, {0, 0}};
    th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr << "Combination of treatments does not work - block 4\n";
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }
    return num_errors;
}

int test_pesticide_temporal_overlap()
{
    int num_errors = 0;
    Scheduler scheduler(Date(2020, 1, 1), Date(2020, 12, 31), StepUnit::Day, 7);
    TestEnvironment environment;
    Raster<double> tr1 = {{1, 1}, {0, 0}};
    Raster<double> tr2 = {{0, 0}, {1, 1}};

    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    Raster<int> zeros(infected.rows(), infected.cols(), 0);
    auto total_hosts = infected + susceptible + resistant;
    std::vector<Raster<int>> exposed;
    std::vector<Raster<int>> mortality_tracker(1, infected);

    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    StandardSingleHostPool host_pool(
        ModelType::SusceptibleInfected,
        bool(mortality_tracker.size()),
        susceptible,
        exposed,
        0,
        infected,
        zeros,
        resistant,
        mortality_tracker,
        zeros,
        total_hosts,
        environment,
        false,
        0,
        false,
        0,
        infected.rows(),
        infected.cols(),
        suitable_cells);

    Treatments<StandardSingleHostPool, Raster<double>> treatments(scheduler);
    treatments.add_treatment(tr1, Date(2020, 5, 1), 30, TreatmentApplication::Ratio);
    treatments.add_treatment(tr2, Date(2020, 5, 20), 30, TreatmentApplication::Ratio);

    unsigned n = scheduler.schedule_action_date(Date(2020, 5, 1));
    treatments.manage(n, host_pool);

    Raster<int> treated = {{0, 0}, {20, 42}};
    Raster<int> inf_treated = {{0, 0}, {16, 40}};
    Raster<int> resist = {{11, 10}, {0, 0}};
    auto th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr
            << "Temporal overlap of pesticide treatments does not work - block 1\n";
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }

    n = scheduler.schedule_action_date(Date(2020, 5, 20));
    treatments.manage(n, host_pool);

    treated = {{0, 0}, {0, 0}};
    inf_treated = {{0, 0}, {0, 0}};
    resist = {{11, 10}, {36, 82}};
    th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr
            << "Temporal overlap of pesticide treatments does not work - block 2\n";
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }

    n = scheduler.schedule_action_date(Date(2020, 6, 1));
    treatments.manage(n, host_pool);

    treated = {{11, 10}, {0, 0}};
    inf_treated = {{0, 0}, {0, 0}};
    resist = {{0, 0}, {36, 82}};
    th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr
            << "Temporal overlap of pesticide treatments does not work - block 3\n";
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }

    n = scheduler.schedule_action_date(Date(2020, 6, 21));
    treatments.manage(n, host_pool);

    treated = {{11, 10}, {36, 82}};
    inf_treated = {{0, 0}, {0, 0}};
    resist = {{0, 0}, {0, 0}};
    th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr
            << "Temporal overlap of pesticide treatments does not work - block 4\n";
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }

    return num_errors;
}

int test_steering()
{
    int num_errors = 0;
    Scheduler scheduler(Date(2020, 1, 1), Date(2020, 12, 31), StepUnit::Day, 7);
    TestEnvironment environment;

    Raster<double> tr1 = {{1, 0.5}, {0.75, 0}};
    Raster<double> tr2 = {{1, 1}, {1, 1}};

    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    Raster<int> zeros(infected.rows(), infected.cols(), 0);
    auto total_hosts = infected + susceptible + resistant;
    std::vector<Raster<int>> exposed;
    std::vector<Raster<int>> mortality_tracker(1, infected);

    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    StandardSingleHostPool host_pool(
        ModelType::SusceptibleInfected,
        bool(mortality_tracker.size()),
        susceptible,
        exposed,
        0,
        infected,
        zeros,
        resistant,
        mortality_tracker,
        zeros,
        total_hosts,
        environment,
        false,
        0,
        false,
        0,
        infected.rows(),
        infected.cols(),
        suitable_cells);
    Treatments<StandardSingleHostPool, Raster<double>> treatments(scheduler);
    treatments.add_treatment(tr1, Date(2020, 5, 1), 0, TreatmentApplication::Ratio);
    treatments.add_treatment(tr2, Date(2020, 6, 1), 7, TreatmentApplication::Ratio);
    unsigned n = scheduler.schedule_action_date(Date(2020, 1, 1));
    treatments.manage(n, host_pool);
    n = scheduler.schedule_action_date(Date(2020, 5, 3));
    treatments.manage(n, host_pool);
    n = scheduler.schedule_action_date(Date(2020, 5, 12));
    treatments.manage(n, host_pool);
    n = scheduler.schedule_action_date(Date(2020, 6, 1));
    treatments.manage(n, host_pool);
    n = scheduler.schedule_action_date(Date(2020, 6, 8));
    treatments.manage(n, host_pool);
    n = scheduler.schedule_action_date(Date(2020, 6, 15));
    treatments.manage(n, host_pool);

    Raster<int> treated = {{0, 5}, {9, 82}};
    Raster<int> inf_treated = {{0, 0}, {0, 0}};
    Raster<int> resist = {{0, 0}, {0, 0}};
    auto th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr << "Steering with treatments does not work - block 1\n";
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }

    susceptible = {{10, 6}, {20, 42}};
    resistant = {{0, 0}, {0, 0}};
    infected = {{1, 4}, {16, 40}};
    mortality_tracker = {infected};
    n = scheduler.schedule_action_date(Date(2020, 1, 1));
    treatments.manage(n, host_pool);
    n = scheduler.schedule_action_date(Date(2020, 5, 3));
    treatments.manage(n, host_pool);
    n = scheduler.schedule_action_date(Date(2020, 5, 12));
    treatments.manage(n, host_pool);
    n = scheduler.schedule_action_date(Date(2020, 6, 1));
    treatments.manage(n, host_pool);
    n = scheduler.schedule_action_date(Date(2020, 6, 8));
    treatments.manage(n, host_pool);
    n = scheduler.schedule_action_date(Date(2020, 6, 15));
    treatments.manage(n, host_pool);

    treated = {{0, 5}, {9, 82}};
    inf_treated = {{0, 0}, {0, 0}};
    resist = {{0, 0}, {0, 0}};
    th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr << "Steering with treatments does not work - block 2\n";
        std::cerr << "Values:\n";
        std::cerr << susceptible << infected << resistant << total_hosts;
        std::cerr << "Differences:\n";
        std::cerr << susceptible - treated << infected - inf_treated
                  << resist - resistant << total_hosts - th_treated;
        num_errors++;
    }

    return num_errors;
}

int test_clear()
{
    int num_errors = 0;
    int num_actions = 0;
    Scheduler scheduler(Date(2020, 1, 1), Date(2020, 12, 31), StepUnit::Day, 7);
    TestEnvironment environment;
    Raster<double> tr1 = {{1, 0.5}, {0.75, 0}};
    Raster<double> tr2 = {{1, 1}, {1, 1}};
    Raster<double> tr3 = {{1, 0}, {1, 0}};

    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    Raster<int> zeros(infected.rows(), infected.cols(), 0);
    auto total_hosts = infected + susceptible + resistant;
    std::vector<Raster<int>> exposed;
    std::vector<Raster<int>> mortality_tracker(1, infected);

    std::vector<std::vector<int>> suitable_cells = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    StandardSingleHostPool host_pool(
        ModelType::SusceptibleInfected,
        bool(mortality_tracker.size()),
        susceptible,
        exposed,
        0,
        infected,
        zeros,
        resistant,
        mortality_tracker,
        zeros,
        total_hosts,
        environment,
        false,
        0,
        false,
        0,
        infected.rows(),
        infected.cols(),
        suitable_cells);
    Treatments<StandardSingleHostPool, Raster<double>> treatments(scheduler);
    treatments.add_treatment(tr1, Date(2020, 5, 1), 0, TreatmentApplication::Ratio);
    treatments.add_treatment(tr2, Date(2020, 6, 1), 7, TreatmentApplication::Ratio);
    treatments.add_treatment(tr3, Date(2020, 6, 8), 7, TreatmentApplication::Ratio);
    unsigned n = scheduler.schedule_action_date(Date(2020, 6, 1));
    treatments.clear_after_step(n);
    n = scheduler.schedule_action_date(Date(2020, 1, 1));
    num_actions += treatments.manage(n, host_pool);
    n = scheduler.schedule_action_date(Date(2020, 5, 3));
    num_actions += treatments.manage(n, host_pool);
    n = scheduler.schedule_action_date(Date(2020, 5, 12));
    num_actions += treatments.manage(n, host_pool);
    n = scheduler.schedule_action_date(Date(2020, 6, 1));
    num_actions += treatments.manage(n, host_pool);
    n = scheduler.schedule_action_date(Date(2020, 6, 8));
    num_actions += treatments.manage(n, host_pool);
    n = scheduler.schedule_action_date(Date(2020, 6, 15));
    num_actions += treatments.manage(n, host_pool);

    Raster<int> treated = {{0, 5}, {9, 82}};
    Raster<int> inf_treated = {{0, 0}, {0, 0}};
    Raster<int> resist = {{0, 0}, {0, 0}};
    auto th_treated = treated + inf_treated + resist;
    if (!(susceptible == treated && infected == inf_treated && resist == resistant
          && total_hosts == th_treated)) {
        std::cerr << "Clearing of treatments does not work - block 1\n";
        std::cerr << susceptible << infected << resistant << total_hosts;
        num_errors++;
    }
    if (num_actions != 3) {
        std::cerr << "Clearing of treatments does not work - block 2\n";
        std::cerr << num_actions << std::endl;
        num_errors++;
    }
    return num_errors;
}

int test_treat_app_from_string()
{
    int num_errors = 0;
    if (treatment_app_enum_from_string("ratio_to_all") != TreatmentApplication::Ratio)
        num_errors++;
    if (treatment_app_enum_from_string("all_infected_in_cell")
        != TreatmentApplication::AllInfectedInCell)
        num_errors++;
    try {
        treatment_app_enum_from_string("invalid_input");
        num_errors++;
    }
    catch (std::invalid_argument&) {
        // OK
    }
    catch (...) {
        num_errors++;
    }
    return num_errors;
}

int main()
{
    int num_errors = 0;

    num_errors += test_application_ratio();
    num_errors += test_application_ratio_without_mortality();
    num_errors += test_application_all_inf();
    num_errors += test_application_ratio_pesticide();
    num_errors += test_application_all_inf_pesticide();
    num_errors += test_combination();
    num_errors += test_pesticide_temporal_overlap();
    num_errors += test_steering();
    num_errors += test_clear();
    num_errors += test_treat_app_from_string();

    (num_errors ? std::cerr : std::cout)
        << "Test treatments number of errors: " << num_errors << "\n";
    return num_errors;
}

#endif  // POPS_TEST
