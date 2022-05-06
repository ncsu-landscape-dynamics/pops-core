
class Hosts
{
    Raster infected;
    vector<Raster> exposed;
    Raster susceptible;
    Raster total_hosts;  // computed on demand
    Raster all_populations;  // computed on demand
    Raster other_individuals;  // maybe environment.other_individuals (non hosts)
    vector<cell> suitable_cells;
    vector<Raster>mortality_tracker_vector;
    Raster died;

    void disperser_to(i, j)
    {
        if (establish) {
            if ("SI") {
                infected(i, j) += 1;
            }
            if ("SEI") {
                exposed.back()(i, j) += 1;
            }
            established_dispersers(i, j) += 1;
        }
    }
    void step_forward() // or infect
    {
        if ("SEI") {
            infected = exposed.front();
            rotate_left_by_one(exposed);
        }
        if (mortality) {
            rotate_left_by_one(mortality_tracker_vector);
        }
    }
    void apply_mortality_at(i, j, values)
    {
        for value in values {
            mortality_tracker_vector[index](i, j) -= value;
            died(i, j) += value;
            if (infected(i, j) > 0)
                infected(i, j) -= value;
            if (total_hosts(i, j) > 0)
                total_hosts(i, j) -= value;
        }
    }
}

class Pests
{
    Raster mobile_dispersers;
    Raster in_cell_dispersers;
    // maybe in soil stored dispersers
}

class Spread
{
    // Generate and kernel part of disperse could be the same,
    // then the disperser moves would be resolved right away,
    // and there would be no need to distinguish between generated dispersers,
    // but it would be much harder to track origin of established dispersers.

    void disperse(Hosts hosts)
    {
        for disperser in dispersers {
            i, j = kernel(disperser);
            hosts.disperser_to(i, j);
        }
    }
    void infect(Hosts hosts)
    {

    }
}

class SurvivalRate
{
    void action(Hosts hosts)
    {
        for i, j in hosts.suitable_cells {
            count = hosts.infected_at(i, j) * survival_rate(i, j);
            hosts.remove_infection_at(i, j, count);
        }
    }
}

class RemoveByTemperature
{
    void action(Hosts hosts)
    {
        for i, j in hosts.suitable_cells {
            if (environment.temperature_at(i, j) < lethal_temperature) {
                count = hosts.infected_at(i, j);
                hosts.remove_infection_at(i, j, count);
            }
        }
    }
}

class PesticideTreatment
{
    void apply_treatment(Hosts hosts)
    {
        for i, j in hosts.suitable_cells {
            susceptible = treatment(i, j) * hosts.get_susceptible_at(i, j);
            infected = get_resistant(treatment(i, j), hosts.infected_at(i, j));
            vector<int> exposed_list;
            for value in hosts.exposed_at(i, j) {
                exposed_list.push_back(get_resistant(treatment(i, j), value));
            }
            hosts.add_resistant_at(i, j, infected, susceptible, exposed_list);  // This could be an apply function with lambda.
        }
    }
    void void end_treatment_effect(Hosts hosts)
    {
        for i, j in hosts.suitable_cells {
            hosts.remove_all_resistant_at(i, j);
        }
    }
}

class Mortality
{
    void action(Hosts host)
    {
        for i, j in hosts.suitable_cells {
            for item in hosts.mortality_at(i, j) {
                if (index == 0)
                    mortality = item;
                else
                    mortality = mortality_rate * item;
                list.push_back(mortality);
            }
            hosts.kill_at(i, j, mortality);
        }
    }
}

class MoveOverpopulatedPests
{
    void action(Hosts hosts)
    {
        // generate_moves as from cell, to cell, count won't work because
        // 1) they need to be airborne first and
        // 2) only count for leaving=infected->susceptible is not tracked (no exposed and resistant)
    }
}

class Model
{
    Hosts hosts;
    hosts.set_infected(...)
    ... // more rasters added here (or all are set in the ctor)

    remove_by_temperature.action(hosts);
    survival_rate.action(hosts);

    move_overpopulated_pests

    // The treatment would be wrapped in Treatments which has manage to decide
    // the time, but it does not need manage_mortality anymore because the
    // loop over mortality tracker is in Hosts.
    pesticide_treatment.apply_treatment(hosts);
    pesticide_treatment.end_treatment_effect(hosts);

    hosts.step_forward();
}
