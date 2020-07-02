/*
 * Tests for the PoPS Config class.
 *
 * Copyright (C) 2020 by the authors.
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

#ifndef POPS_CONFIG_HPP
#define POPS_CONFIG_HPP

#include "scheduling.hpp"

#include <vector>

namespace pops {

class Config
{
public:
    Config()
        :
          generate_stochasticity(true),
          establishment_stochasticity(true),
          use_lethal_temperature(false),
          use_anthropogenic_kernel(false),
          use_treatments(false),
          use_mortality(false)
    {}
    // Seed
    int random_seed;
    // Size
    int rows;
    int cols;
    double ew_res;
    double ns_res;
    // TODO: Time
    //int steps;
    // Reduced stochasticity
    bool generate_stochasticity;
    bool establishment_stochasticity;
    double establishment_probability;
    // Temperature
    bool use_lethal_temperature;
    double lethal_temperature;
    bool weather;
    double reproductive_rate;
    // SI/SEI
    std::string model_type;
    int latency_period_steps;
    // Kernels
    std::string natural_kernel_type;
    double natural_scale;
    std::string natural_direction;
    double natural_kappa;
    bool use_anthropogenic_kernel;
    double percent_natural_dispersal;
    std::string anthro_kernel_type;
    double anthro_scale;
    std::string anthro_direction;
    double anthro_kappa;
    // Treatments
    bool use_treatments;
    // Mortality
    bool use_mortality;
    double mortality_rate;
    int first_mortality_year;  // TODO: document that it starts at 1, not 0

    std::string step_unit;
    unsigned step_num_units;

    int season_start_month{0};
    int season_end_month{0};

    std::string output_frequency;
    unsigned output_frequency_n;
    int lethal_temperature_month;

    void create_schedules()
    {
        scheduler_ = Scheduler(Date(date_start_), Date(date_end_), step_unit_enum_from_string(step_unit), step_num_units);
        spread_schedule_ = scheduler_.schedule_spread(Season(season_start_month, season_end_month));
        output_schedule_ = output_schedule_from_string(scheduler_, output_frequency, output_frequency_n);
        mortality_schedule_ = scheduler_.schedule_action_end_of_year();
        lethal_schedule_ = scheduler_.schedule_action_yearly(lethal_temperature_month, 1);
        spread_rate_schedule_ = scheduler_.schedule_action_end_of_year();
        schedules_created_ = true;
    }

    const Scheduler& scheduler() const
    {
        if (!schedules_created_)
            throw std::logic_error("Schedules were not created before calling scheduler()");
        return scheduler_;
    }

    const std::vector<bool>& spread_schedule() const
    {
        if (!schedules_created_)
            throw std::logic_error("Schedules were not created before calling spread_schedule()");
        return spread_schedule_;
    }

    const std::vector<bool>& mortality_schedule() const
    {
        if (!schedules_created_)
            throw std::logic_error("Schedules were not created before calling mortality_schedule()");
        return mortality_schedule_;
    }

    const std::vector<bool>& lethal_schedule() const
    {
        if (!schedules_created_)
            throw std::logic_error("Schedules were not created before calling lethal_schedule()");
        return lethal_schedule_;
    }

    const std::vector<bool>& spread_rate_schedule() const
    {
        if (!schedules_created_)
            throw std::logic_error("Schedules were not created before calling spread_rate_schedule()");
        return spread_rate_schedule_;
    }

    const std::vector<bool>& output_schedule() const
    {
        if (!schedules_created_)
            throw std::logic_error("Schedules were not created before calling output_schedule()");
        return output_schedule_;
    }

    unsigned num_mortality_years()
    {
        if (!schedules_created_)
            throw std::logic_error("Schedules were not created before calling num_mortality_years()");
        return get_number_of_scheduled_actions(mortality_schedule_);
    }

    unsigned num_lethal()
    {
        if (!schedules_created_)
            throw std::logic_error("Schedules were not created before calling num_lethal()");
        get_number_of_scheduled_actions(lethal_schedule_);
    }

    unsigned rate_num_years()
    {
        if (!schedules_created_)
            throw std::logic_error("Schedules were not created before calling rate_num_years()");
        get_number_of_scheduled_actions(spread_rate_schedule_);
    }

    const std::string& date_end() const
    {
        return date_end_;
    }

    void set_date_end(const std::string &date_end)
    {
        date_end_ = date_end;
    }

    const std::string& date_start() const
    {
        return date_start_;
    }

    void set_date_start(const std::string &date_start)
    {
        date_start_ = date_start;
    }

private:
    std::string date_start_;
    std::string date_end_;

    Scheduler scheduler_{Date("0-01-01"), Date("0-01-02"), StepUnit::Day, 1};
    bool schedules_created_{false};

    std::vector<bool> spread_schedule_;
    std::vector<bool> output_schedule_;
    std::vector<bool> mortality_schedule_;
    std::vector<bool> lethal_schedule_;
    std::vector<bool> spread_rate_schedule_;
};

}  // namespace pops

#endif  // POPS_CONFIG_HPP
