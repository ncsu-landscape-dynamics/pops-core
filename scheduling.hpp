/*
 * SOD model - scheduling simulation steps
 *
 * Copyright (C) 2015-2019 by the authors.
 *
 * Authors: Anna Petrasova
 *          Vaclav Petras
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */


#ifndef POPS_SCHEDULING_HPP
#define POPS_SCHEDULING_HPP

#include <iostream>
#include <vector>

#include "date.hpp"

namespace pops {

/*! Representation and manipulation of a date for the simulation.
 *
 */

/*!
 * 
 */
class Step
{
public:
    Step(Date start_date, Date end_date)
        : start_date_(start_date), end_date_(end_date)
    {}
    inline Date start_date(){return start_date_;}
    inline Date end_date(){return end_date_;}
    inline friend std::ostream& operator<<(std::ostream& os, const Step &step);
private:
    Date start_date_;
    Date end_date_;
};
std::ostream& operator<<(std::ostream& os, const Step &step)
{
    os << step.start_date_ << " - " << step.end_date_;
    return os;
}

/**
 * @brief Enum for step unit
 */
enum class StepUnit {
    Day, Month
};

class Scheduler
{
public:
    Scheduler(const Date start, const Date end, StepUnit spread_unit, unsigned spread_num_units)
        :
          start_(start),
          end_(end),
          spread_unit_(spread_unit),
          spread_num_units_(spread_num_units)
    {
        if (start >= end)
            throw std::invalid_argument("Start date must be before end date");
        Date d(start);
        increase_date(d);
        if (d > end)
            throw std::invalid_argument("There must be at least one step between start and end date");
        if (spread_unit == StepUnit::Month && start.day() != 1)
            throw std::invalid_argument("If step unit is month, start date must start the first day of a month");
        
        Date date(start_);
        unsigned step = 0;
        while (date < end_) {
            Date start(date);
            increase_date(date);
            Date end(date);
            end.subtract_day();
            steps.push_back(Step(start, end));
            step++;
        }
        num_steps = step;
    }

    unsigned get_num_steps() {
        return num_steps;
    }

    std::vector<bool> schedule_spread(const Season &season) {
        std::vector<bool> schedule;
        schedule.reserve(num_steps);
        for (Step step : steps) {
            if (season.month_in_season(step.start_date().month()) || season.month_in_season(step.end_date().month()))
                schedule.push_back(true);
            else
                schedule.push_back(false);
        }
        return schedule;
    }
    std::vector<bool> schedule_action_yearly(int month, int day) {
        std::vector<bool> schedule;
        schedule.reserve(num_steps);
        for (Step step : steps) {
            Date st = step.start_date();
            Date end = step.end_date();
            /* this doesn't handle cases where
               there is change in year within interval,
               but that doesn't happen with current implementation */
            Date test(st.year(), month, day);
            if ((test >= st && test <= end))
                schedule.push_back(true);
            else
                schedule.push_back(false);
        }
        return schedule;
    }

    std::vector<bool> schedule_action_end_of_year() {
        std::vector<bool> schedule;
        schedule.reserve(num_steps);
        for (Step step : steps) {
            if (step.end_date().is_last_day_of_year())
                schedule.push_back(true);
            else
                schedule.push_back(false);
        }
        return schedule;
    }

    std::vector<bool> schedule_action_nsteps(unsigned n_steps) {
        std::vector<bool> schedule;
        schedule.reserve(num_steps);
        for (unsigned i = 0; i < num_steps; i++) {
            if ((i + 1) % n_steps == 0)
                schedule.push_back(true);
            else
                schedule.push_back(false);
        }
        return schedule;
    }

    std::vector<bool> schedule_action_monthly() {
        std::vector<bool> schedule;
        schedule.reserve(num_steps);
        for (Step step : steps) {
            Date st = step.start_date();
            Date end = step.end_date();
            if (st.month() != end.month() || end.is_last_day_of_month())
                schedule.push_back(true);
            else
                schedule.push_back(false);
        }
        return schedule;
    }

    unsigned schedule_action_date(const Date &date) {
        for (unsigned i = 0; i < num_steps; i++) {
            if (date >= steps[i].start_date() && date <= steps[i].end_date())
                return i;
        }
        throw std::invalid_argument("Date is outside of schedule");
    }
    void debug_schedule(std::vector<bool> &schedule) {
        for (unsigned i = 0; i < num_steps; i++)
            std::cout << steps[i] << ": " << (schedule.at(i) ? "true" : "false") << std::endl;
    }
    void debug_schedule(unsigned n) {
        for (unsigned i = 0; i < num_steps; i++)
            std::cout << steps[i] << ": " << (n == i ? "true" : "false") << std::endl;
    }

private:
    Date start_;
    Date end_;
    StepUnit spread_unit_;
    unsigned spread_num_units_;
    std::vector<Step> steps;
    unsigned num_steps;

    void increase_date(Date &date) {
        if (spread_unit_ == StepUnit::Day) {
            date.increased_by_days(spread_num_units_);
        }
        else if (spread_unit_ == StepUnit::Month) {
            for (unsigned i = 0; i < spread_num_units_; i++) {
                date.increased_by_month();
            }
        }
    }
};



} // namespace pops

#endif // POPS_SCHEDULING_HPP
