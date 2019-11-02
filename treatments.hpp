/*
 * PoPS model - treatments
 *
 * Copyright (C) 2015-2019 by the authors.
 *
 * Authors: Anna Petrasova <akratoc gmail com>
 *          Vaclav Petras <wenzeslaus gmail com>
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_TREATMENTS_HPP
#define POPS_TREATMENTS_HPP

#include "raster.hpp"
#include "date.hpp"

#include <map>
#include <vector>

namespace pops {


/**
 * @brief The enum to decide how treatment is applied
 */
enum class TreatmentApplication {
    Ratio,  ///< A ratio is applied to all treated rasters
    AllInfectedInCell  ///< All infected individuals are removed, rest by ratio
};


template<typename IntegerRaster, typename FloatRaster>
class AbstractTreatment
{
public:
    virtual Date get_start() = 0;
    virtual Date get_end() = 0;
    virtual void reset(const Date& date) = 0;
    virtual bool should_start(const Date& date) = 0;
    virtual bool should_end(const Date& date) = 0;
    virtual bool should_apply_mortality() = 0;
    virtual void apply_treatment(IntegerRaster& infected, IntegerRaster& susceptible, IntegerRaster& resistant) = 0;
    virtual void end_treatment(IntegerRaster& susceptible, IntegerRaster& resistant) = 0;
    virtual void apply_treatment_mortality(IntegerRaster& infected) = 0;
    virtual ~AbstractTreatment() {}
};


template<typename IntegerRaster, typename FloatRaster>
class BaseTreatment : public AbstractTreatment<IntegerRaster, FloatRaster>
{
protected:
    Date start;
    Date end;
    bool active;
    bool apply_mortality;
    FloatRaster map;
    TreatmentApplication application;
public:
    BaseTreatment(const FloatRaster& map, const Date& start,
                  TreatmentApplication treatment_application):
        start(start), end(start), active(false), apply_mortality(false), map(map), application(treatment_application)
    {}
    Date get_start() {return start;}
    Date get_end() {return end;}
    bool should_apply_mortality() override
    {
        return apply_mortality;
    }
    void apply_treatment_mortality(IntegerRaster& infected) override
    {
        for(unsigned i = 0; i < infected.rows(); i++)
            for(unsigned j = 0; j < infected.cols(); j++) {
                if (application == TreatmentApplication::Ratio) {
                    infected(i, j) = infected(i, j) * map(i, j);
                }
                else if (application == TreatmentApplication::AllInfectedInCell) {
                    infected(i, j) = map(i, j) ? 0 : infected(i, j);
                }
            }
        apply_mortality = false;
    }
};


template<typename IntegerRaster, typename FloatRaster>
class SimpleTreatment : public BaseTreatment<IntegerRaster, FloatRaster>
{
public:
    SimpleTreatment(const FloatRaster& map, const Date& start,
                    TreatmentApplication treatment_application):
        BaseTreatment<IntegerRaster, FloatRaster>(map, start, treatment_application)
    {}
    bool should_start(const Date& date) override
    {
        if (!this->active && date >= this->start)
            return true;
        return false;
    }
    bool should_end(const Date&) override
    {
        return false;
    }
    void reset(const Date& date) override {
        if (this->get_start() >= date)
            this->active = false;
    }
    void apply_treatment(IntegerRaster& infected, IntegerRaster& susceptible, IntegerRaster& ) override
    {
        for(unsigned i = 0; i < infected.rows(); i++)
            for(unsigned j = 0; j < infected.cols(); j++) {
                if (this->application == TreatmentApplication::Ratio) {
                    infected(i, j) = infected(i, j) - (infected(i, j) * this->map(i, j));
                }
                else if (this->application == TreatmentApplication::AllInfectedInCell) {
                    infected(i, j) = this->map(i, j) ? 0 : infected(i, j);
                }
                susceptible(i, j) = susceptible(i, j) - (susceptible(i, j) * this->map(i, j));
            }
        this->active = true;
    }
    void end_treatment(IntegerRaster&, IntegerRaster&) override
    {
        return;
    }
};


template<typename IntegerRaster, typename FloatRaster>
class PesticideTreatment : public BaseTreatment<IntegerRaster, FloatRaster>
{
public:
    PesticideTreatment(const FloatRaster& map, const Date& start, int num_days,
                       TreatmentApplication treatment_application):
        BaseTreatment<IntegerRaster, FloatRaster>(map, start, treatment_application)
    {
        this->end.increased_by_days(num_days);
    }
    bool should_start(const Date& date) override
    {
        if (!this->active && date >= this->start && date < this->end)
            return true;
        return false;
    }
    bool should_end(const Date& date) override
    {
        if (this->active && date >= this->end)
            return true;
        return false;
    }
    void reset(const Date&) override {
        this->active = false;
    }
    void apply_treatment(IntegerRaster& infected, IntegerRaster& susceptible, IntegerRaster& resistant) override
    {
        for(unsigned i = 0; i < infected.rows(); i++)
            for(unsigned j = 0; j < infected.cols(); j++) {
                int infected_resistant;
                int susceptible_resistant = susceptible(i, j) * this->map(i, j);
                if (this->application == TreatmentApplication::Ratio) {
                    infected_resistant = infected(i, j) * this->map(i, j);
                }
                else if (this->application == TreatmentApplication::AllInfectedInCell) {
                    infected_resistant = this->map(i, j) ? infected(i, j): 0;
                }
                infected(i, j) -= infected_resistant;
                resistant(i, j) = infected_resistant + susceptible_resistant;
                susceptible(i, j) -= susceptible_resistant;
            }
        this->active = true;
        this->apply_mortality = true;
    }
    void end_treatment(IntegerRaster& susceptible, IntegerRaster& resistant) override
    {
        for(unsigned i = 0; i < resistant.rows(); i++)
            for(unsigned j = 0; j < resistant.cols(); j++) {
                if (this->map(i, j) > 0){
                    susceptible(i, j) += resistant(i, j);
                    resistant(i, j) = 0;
                }
            }
        this->active = false;
    }
};


template<typename IntegerRaster, typename FloatRaster>
class Treatments
{
private:
    std::vector<AbstractTreatment<IntegerRaster, FloatRaster>*> treatments;
public:
    ~Treatments()
    {
        for (auto item : treatments)
        {
            delete item;
        }
    }
    void add_treatment(const FloatRaster& map, const Date& start_date, int num_days, TreatmentApplication treatment_application)
    {
        if (num_days == 0)
            treatments.push_back(new SimpleTreatment<IntegerRaster, FloatRaster>(map, start_date, treatment_application));
        else
            treatments.push_back(new PesticideTreatment<IntegerRaster, FloatRaster>(map, start_date, num_days, treatment_application));
    }
    bool manage(const Date& current, IntegerRaster& infected,
                IntegerRaster& susceptible, IntegerRaster& resistant)
    {
        bool applied = false;
        for (unsigned i = 0; i < treatments.size(); i++) {
            if (treatments[i]->should_start(current)) {
                treatments[i]->apply_treatment(infected, susceptible, resistant);
                applied = true;
            }
            else if (treatments[i]->should_end(current)) {
                treatments[i]->end_treatment(susceptible, resistant);
                applied = true;
            }
        }
        return applied;
    }
    bool manage_mortality(IntegerRaster& infected)
    {
        bool applied = false;
        for (unsigned i = 0; i < treatments.size(); i++)
            if (treatments[i]->should_apply_mortality()) {
                treatments[i]->apply_treatment_mortality(infected);
                applied = true;
            }
        return applied;
    }
    /* this needs to be called when going back
       when doing computational steering */
    void reset(const Date& date)
    {
        for (auto item : treatments)
        {
            item->reset(date);
            item->get_start();
        }
    }
    void clear_after_date(const Date& date)
    {
        for(auto& treatment : treatments)
        {
            if (treatment->get_start() > date)
            {
                delete treatment;
                treatment = nullptr;
            }
        }
        treatments.erase(std::remove(treatments.begin(), treatments.end(), nullptr), treatments.end());
    }
};

}
#endif // POPS_TREATMENTS_HPP

