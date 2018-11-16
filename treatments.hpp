/*
 * PoPS model - treatments
 *
 * Copyright (C) 2015-2018 by the authors.
 *
 * Authors: Anna Petrasova
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

#include <map>

namespace pops {

template<typename IntegerRaster>
class Treatments
{
private:
    std::map<int, IntegerRaster> treatments;
public:
    void add_treatment(int year, const IntegerRaster &map)
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
    void apply_treatment(int year, IntegerRaster &host)
    {
        if (treatments.find(year) != treatments.end()) {
            host = host - (host * treatments[year]);
        }
        // otherwise no treatment for that year
    }
};

}
#endif // POPS_TREATMENTS_HPP

