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

#include "raster.hpp"
#include "treatments.hpp"


using namespace pops;


int test_application_ratio()
{
    int num_errors = 0;
    Treatments<Raster<int>, Raster<double>> treatments;
    Raster<double> tr1 = {{1, 0.5}, {0.75, 0}};
    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    treatments.add_treatment(tr1, Date(2020, 1, 1), 0, TreatmentApplication::Ratio);
    treatments.manage(Date(2020, 1, 1), infected, susceptible, resistant);

    Raster<int> treated = {{0, 3}, {5, 42}};
    Raster<int> inf_treated = {{0, 2}, {4, 40}};
    if (!(susceptible == treated && infected == inf_treated)) {
        std::cout << "Treatment with ratio app does not work" << std::endl;
        std::cout << susceptible << infected;
        num_errors++;
    }
    return num_errors;
}

int test_application_all_inf()
{
    int num_errors = 0;
    Treatments<Raster<int>, Raster<double>> treatments;
    Raster<double> tr1 = {{1, 0.5}, {0.75, 0}};
    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    treatments.add_treatment(tr1, Date(2020, 1, 1), 0, TreatmentApplication::AllInfectedInCell);
    treatments.manage(Date(2020, 1, 1), infected, susceptible, resistant);

    Raster<int> treated = {{0, 3}, {5, 42}};
    Raster<int> inf_treated = {{0, 0}, {0, 40}};
    if (!(susceptible == treated && infected == inf_treated)) {
        std::cout << "Treatment with AllInfectedInCell app does not work" << std::endl;
        std::cout << susceptible << infected;
        num_errors++;
    }
    return num_errors;
}

int test_application_ratio_pesticide()
{
    int num_errors = 0;
    Treatments<Raster<int>, Raster<double>> treatments;
    Raster<double> tr1 = {{1, 0.5}, {0.75, 0}};
    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    treatments.add_treatment(tr1, Date(2020, 5, 1), 7, TreatmentApplication::Ratio);
    treatments.manage(Date(2020, 1, 1), infected, susceptible, resistant);

    Raster<int> treated = {{10, 6}, {20, 42}};
    Raster<int> inf_treated = {{1, 4}, {16, 40}};
    Raster<int> resist = {{0, 0}, {0, 0}};
    if (!(susceptible == treated && infected == inf_treated && resist == resistant)) {
        std::cout << "Treatment with pesticide and ratio app does not work" << std::endl;
        std::cout << susceptible << infected << resistant;
        num_errors++;
    }
    treatments.manage(Date(2020, 5, 3), infected, susceptible, resistant);

    treated = {{0, 3}, {5, 42}};
    inf_treated = {{0, 2}, {4, 40}};
    resist = {{11, 5}, {27, 0}};
    if (!(susceptible == treated && infected == inf_treated && resist == resistant)) {
        std::cout << "Treatment with pesticide and ratio app does not work" << std::endl;
        std::cout << susceptible << infected << resistant;
        num_errors++;
    }
    treatments.manage(Date(2020, 5, 8), infected, susceptible, resistant);

    treated = {{11, 8}, {32, 42}};
    resist = {{0, 0}, {0, 0}};
    if (!(susceptible == treated && infected == inf_treated && resist == resistant)) {
        std::cout << "Treatment with pesticide and ratio app does not work" << std::endl;
        std::cout << susceptible << infected << resistant;
        num_errors++;
    }
    return num_errors;
}


int test_application_all_inf_pesticide()
{
    int num_errors = 0;
    Treatments<Raster<int>, Raster<double>> treatments;
    Raster<double> tr1 = {{1, 0.5}, {0.75, 0}};
    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    treatments.add_treatment(tr1, Date(2020, 5, 1), 7, TreatmentApplication::AllInfectedInCell);
    treatments.manage(Date(2020, 1, 1), infected, susceptible, resistant);

    Raster<int> treated = {{10, 6}, {20, 42}};
    Raster<int> inf_treated = {{1, 4}, {16, 40}};
    Raster<int> resist = {{0, 0}, {0, 0}};
    if (!(susceptible == treated && infected == inf_treated && resist == resistant)) {
        std::cout << "Treatment with pesticide and AllInfectedInCell app does not work" << std::endl;
        std::cout << susceptible << infected << resistant;
        num_errors++;
    }
    treatments.manage(Date(2020, 5, 3), infected, susceptible, resistant);

    treated = {{0, 3}, {5, 42}};
    inf_treated = {{0, 0}, {0, 40}};
    resist = {{11, 7}, {31, 0}};
    if (!(susceptible == treated && infected == inf_treated && resist == resistant)) {
        std::cout << "Treatment with pesticide and AllInfectedInCell app does not work" << std::endl;
        std::cout << susceptible << infected << resistant;
        num_errors++;
    }
    treatments.manage(Date(2020, 5, 8), infected, susceptible, resistant);

    treated = {{11, 10}, {36, 42}};
    resist = {{0, 0}, {0, 0}};
    if (!(susceptible == treated && infected == inf_treated && resist == resistant)) {
        std::cout << "Treatment with pesticide and AllInfectedInCell app does not work" << std::endl;
        std::cout << susceptible << infected << resistant;
        num_errors++;
    }
    return num_errors;
}


int test_combination()
{
    int num_errors = 0;
    Treatments<Raster<int>, Raster<double>> treatments;
    Raster<double> tr1 = {{1, 0.5}, {0.75, 0}};
    Raster<double> tr2 = {{1, 1}, {1, 1}};

    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    treatments.add_treatment(tr1, Date(2020, 5, 1), 0, TreatmentApplication::Ratio);
    treatments.add_treatment(tr2, Date(2020, 6, 1), 7, TreatmentApplication::Ratio);
    treatments.manage(Date(2020, 1, 1), infected, susceptible, resistant);

    Raster<int> treated = {{10, 6}, {20, 42}};
    Raster<int> inf_treated = {{1, 4}, {16, 40}};
    Raster<int> resist = {{0, 0}, {0, 0}};
    if (!(susceptible == treated && infected == inf_treated && resist == resistant)) {
        std::cout << "Combination of treatments does not work" << std::endl;
        std::cout << susceptible << infected << resistant;
        num_errors++;
    }
    treatments.manage(Date(2020, 5, 3), infected, susceptible, resistant);

    treated = {{0, 3}, {5, 42}};
    inf_treated = {{0, 2}, {4, 40}};
    resist = {{0, 0}, {0, 0}};
    if (!(susceptible == treated && infected == inf_treated && resist == resistant)) {
        std::cout << "Combination of treatments does not work" << std::endl;
        std::cout << susceptible << infected << resistant;
        num_errors++;
    }
    treatments.manage(Date(2020, 6, 3), infected, susceptible, resistant);

    treated = {{0, 0}, {0, 0}};
    inf_treated = {{0, 0}, {0, 0}};
    resist = {{0, 5}, {9, 82}};
    if (!(susceptible == treated && infected == inf_treated && resist == resistant)) {
        std::cout << "Combination of treatments does not work" << std::endl;
        std::cout << susceptible << infected << resistant;
        num_errors++;
    }
    treatments.manage(Date(2020, 6, 8), infected, susceptible, resistant);

    treated = {{0, 5}, {9, 82}};
    inf_treated = {{0, 0}, {0, 0}};
    resist = {{0, 0}, {0, 0}};
    if (!(susceptible == treated && infected == inf_treated && resist == resistant)) {
        std::cout << "Combination of treatments does not work" << std::endl;
        std::cout << susceptible << infected << resistant;
        num_errors++;
    }
    return num_errors;
}

int test_steering()
{
    int num_errors = 0;
    Treatments<Raster<int>, Raster<double>> treatments;
    Raster<double> tr1 = {{1, 0.5}, {0.75, 0}};
    Raster<double> tr2 = {{1, 1}, {1, 1}};

    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    treatments.add_treatment(tr1, Date(2020, 5, 1), 0, TreatmentApplication::Ratio);
    treatments.add_treatment(tr2, Date(2020, 6, 1), 7, TreatmentApplication::Ratio);
    treatments.manage(Date(2020, 1, 1), infected, susceptible, resistant);
    treatments.manage(Date(2020, 5, 3), infected, susceptible, resistant);
    treatments.manage(Date(2020, 5, 12), infected, susceptible, resistant);
    treatments.manage(Date(2020, 6, 1), infected, susceptible, resistant);
    treatments.manage(Date(2020, 6, 3), infected, susceptible, resistant);
    treatments.manage(Date(2020, 6, 8), infected, susceptible, resistant);
    treatments.manage(Date(2020, 6, 10), infected, susceptible, resistant);

    Raster<int> treated = {{0, 5}, {9, 82}};
    Raster<int> inf_treated = {{0, 0}, {0, 0}};
    Raster<int> resist = {{0, 0}, {0, 0}};
    if (!(susceptible == treated && infected == inf_treated && resist == resistant)) {
        std::cout << "Combination of treatments does not work" << std::endl;
        std::cout << susceptible << infected << resistant;
        num_errors++;
    }

    susceptible = {{10, 6}, {20, 42}};
    resistant = {{0, 0}, {0, 0}};
    infected = {{1, 4}, {16, 40}};
    treatments.reset(Date(2020, 1, 1));
    treatments.manage(Date(2020, 1, 1), infected, susceptible, resistant);
    treatments.manage(Date(2020, 5, 3), infected, susceptible, resistant);
    treatments.manage(Date(2020, 5, 12), infected, susceptible, resistant);
    treatments.manage(Date(2020, 6, 1), infected, susceptible, resistant);
    treatments.manage(Date(2020, 6, 3), infected, susceptible, resistant);
    treatments.manage(Date(2020, 6, 8), infected, susceptible, resistant);
    treatments.manage(Date(2020, 6, 10), infected, susceptible, resistant);

    treated = {{0, 5}, {9, 82}};
    inf_treated = {{0, 0}, {0, 0}};
    resist = {{0, 0}, {0, 0}};
    if (!(susceptible == treated && infected == inf_treated && resist == resistant)) {
        std::cout << "Combination of treatments does not work" << std::endl;
        std::cout << susceptible << infected << resistant;
        num_errors++;
    }

    return num_errors;
}

int test_clear()
{
    int num_errors = 0;
    int num_actions = 0;
    Treatments<Raster<int>, Raster<double>> treatments;
    Raster<double> tr1 = {{1, 0.5}, {0.75, 0}};
    Raster<double> tr2 = {{1, 1}, {1, 1}};
    Raster<double> tr3 = {{1, 0}, {1, 0}};

    Raster<int> susceptible = {{10, 6}, {20, 42}};
    Raster<int> resistant = {{0, 0}, {0, 0}};
    Raster<int> infected = {{1, 4}, {16, 40}};
    treatments.add_treatment(tr1, Date(2020, 5, 1), 0, TreatmentApplication::Ratio);
    treatments.add_treatment(tr2, Date(2020, 6, 1), 7, TreatmentApplication::Ratio);
    treatments.add_treatment(tr2, Date(2020, 6, 8), 7, TreatmentApplication::Ratio);
    treatments.clear_after_date(Date(2020, 6, 1));
    num_actions += treatments.manage(Date(2020, 1, 1), infected, susceptible, resistant);
    num_actions += treatments.manage(Date(2020, 5, 3), infected, susceptible, resistant);
    num_actions += treatments.manage(Date(2020, 5, 12), infected, susceptible, resistant);
    num_actions += treatments.manage(Date(2020, 6, 1), infected, susceptible, resistant);
    num_actions += treatments.manage(Date(2020, 6, 3), infected, susceptible, resistant);
    num_actions += treatments.manage(Date(2020, 6, 8), infected, susceptible, resistant);
    num_actions += treatments.manage(Date(2020, 6, 10), infected, susceptible, resistant);

    Raster<int> treated = {{0, 5}, {9, 82}};
    Raster<int> inf_treated = {{0, 0}, {0, 0}};
    Raster<int> resist = {{0, 0}, {0, 0}};
    if (!(susceptible == treated && infected == inf_treated && resist == resistant)) {
        std::cout << "Clearing of treatments does not work" << std::endl;
        std::cout << susceptible << infected << resistant;
        num_errors++;
    }
    if (num_actions != 3) {
        std::cout << "Clearing of treatments does not work" << std::endl;
        std::cout << num_actions << std::endl;
        num_errors++;
    }
    return num_errors;
}

int main()
{
    int num_errors = 0;

    num_errors += test_application_ratio();
    num_errors += test_application_all_inf();
    num_errors += test_application_ratio_pesticide();
    num_errors += test_application_all_inf_pesticide();
    num_errors += test_combination();
    num_errors += test_steering();
    num_errors += test_clear();

    return num_errors;
}

#endif  // POPS_TEST
