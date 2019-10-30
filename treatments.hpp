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
 * @brief The enum to decide how threatment is applied
 */
enum class TreatmentApplication {
    Ratio,  ///< A ratio is applied to all treated rasters
    AllInfectedInCell  ///< All infected individuals are removed, rest by ratio
};

template<typename IntegerRaster, typename FloatRaster>
class Treatments
{
private:
    std::map<int, FloatRaster> treatments;
    TreatmentApplication application;
public:
    Treatments(
        TreatmentApplication treatment_application = TreatmentApplication::Ratio) :
        application(treatment_application)
    {}
    void add_treatment(int year, const FloatRaster &map)
    {
        treatments[year] = map;
    }
    void clear_all()
    {
        treatments.erase(treatments.begin(), treatments.end());
    }
    void clear_after_year(int year)
    {
        for (auto it = treatments.begin(); it != treatments.end();) {
            if (it->first > year) {
                treatments.erase(it++);
            }
            else {
                ++it;
            }
        }
    }
    void apply_treatment_host(int year, IntegerRaster &infected, IntegerRaster &susceptible)
    {
        // this expression fails in rcpp
        // host = host - (host * treatments[year]);
        if (treatments.find(year) != treatments.end()) {
            for(unsigned i = 0; i < infected.rows(); i++)
                for(unsigned j = 0; j < infected.cols(); j++) {
                    if (application == TreatmentApplication::Ratio) {
                        infected(i, j) = infected(i, j) - (infected(i, j) * treatments[year](i, j));
                        susceptible(i, j) = susceptible(i, j) - (susceptible(i, j) * treatments[year](i, j));
                    }
                    else if (application == TreatmentApplication::AllInfectedInCell) {
                        infected(i, j) = treatments[year](i, j) ? 0 : infected(i, j);
                        susceptible(i, j) = susceptible(i, j) - (susceptible(i, j) * treatments[year](i, j));
                    }
                }
        }
        // otherwise no treatment for that year
    }
    void apply_treatment_infected(int year, IntegerRaster &infected)
    {
        if (treatments.find(year) != treatments.end()) {
            for(unsigned i = 0; i < infected.rows(); i++)
                for(unsigned j = 0; j < infected.cols(); j++)
                    if (application == TreatmentApplication::Ratio) {
                        infected(i, j) = infected(i, j) - (infected(i, j) * treatments[year](i, j));
                    }
                    else if (application == TreatmentApplication::AllInfectedInCell) {
                        infected(i, j) = treatments[year](i, j) ? 0 : infected(i, j);
                    }
        }
    }
};

template<typename IntegerRaster, typename FloatRaster>
class PesticideTreatment
{
private:
    Date start;
    Date end;
    bool resistant;
    bool apply_mortality;
    FloatRaster map;
    TreatmentApplication application;
public:
    Treatment(const FloatRaster &map, const Date &start, int num_days,
              TreatmentApplication treatment_application):
        map(map), start(start), resistant(false), apply_mortality(false), application(treatment_application)
    {
        Date end_date = Date(start);
        end_date.increased_by_days(num_days);
    }
    Treatment() = delete;

    bool should_start(const Date &date)
    {
        if (!resistant && date >= start && date < end)
            return true;
        return false;
    }
    bool should_end(const Date &date)
    {
        if (resistant && date >= end)
            return true;
        return false;
    }
    bool should_apply_mortality()
    {
        return apply_mortality;
    }

    void reset_resistance(){
        resistant = false;
    }

    void apply_treatment(IntegerRaster &infected, IntegerRaster &susceptible, IntegerRaster &resistant)
    {
        for(unsigned i = 0; i < infected.rows(); i++)
            for(unsigned j = 0; j < infected.cols(); j++) {
                int infected_resistant;
                int susceptible_resistant = susceptible(i, j) * map(i, j);
                if (application == TreatmentApplication::Ratio) {
                    infected_resistant = infected(i, j) * map(i, j);
                }
                else if (application == TreatmentApplication::AllInfectedInCell) {
                    infected_resistant = map(i, j) ? infected(i, j): 0;
                }
                infected(i, j) -= infected_resistant;
                resistant(i, j) = infected_resistant + susceptible_resistant;
                susceptible(i, j) -= susceptible_resistant;
            }
        resistant = true;
        apply_mortality = true;
    }
    void end_treatment(IntegerRaster &susceptible, IntegerRaster &resistant)
    {
        for(unsigned i = 0; i < resistant.rows(); i++)
            for(unsigned j = 0; j < resistant.cols(); j++) {
                if (map(i, j) > 0){
                    susceptible(i, j) += resistant(i, j);
                    resistant(i, j) = 0;
                }
            }
        resistant = false;
    }
    void apply_treatment_mortality(IntegerRaster &infected)
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
class PesticideTreatments
{
private:
    std::vector<PesticideTreatment<IntegerRaster, FloatRaster>> treatments;
    TreatmentApplication application;

public:
    PesticideTreatments(TreatmentApplication treatment_application = TreatmentApplication::Ratio) :
        application(treatment_application)
    {}
    void add_treatment(const FloatRaster &map, const Date &start_date, int num_days)
    {
        treatments.emplace_back(map, start_date, num_days, application);
    }
    /* this needs to be called when going back
       when doing computational steering */
    void reset_resistance()
    {
        for (unsigned i = 0; i < treatments.size(); i++)
            treatments[i].reset_resistance();
    }

    bool pesticide_management(const Date &current, IntegerRaster &infected,
                              IntegerRaster &susceptible, IntegerRaster &resistant)
    {
        bool applied = false;
        for (unsigned i = 0; i < treatments.size(); i++) {
            if (treatments[i].should_start(current)) {
                treatments[i].apply_treatment(infected, susceptible, resistant);
                applied = true;
            }
            else if (treatments[i].should_end(current)) {
                treatments[i].end_treatment(susceptible, resistant);
                applied = true;
            }
        }
        return applied;
    }

    bool pesticide_management_mortality(IntegerRaster &infected)
    {
        bool applied = false;
        for (unsigned i = 0; i < treatments.size(); i++)
            if (treatments[i].should_apply_mortality()) {
                treatments[i].apply_treatment_mortality(infected);
                applied = true;
            }
        return applied;
    }
};

};
#endif // POPS_TREATMENTS_HPP

