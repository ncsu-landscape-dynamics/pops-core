#include <vector>

namespace pops {

template<typename IntegerRaster>
class Host
{
    IntegerRaster infected;
    std::vector<IntegerRaster> exposed;
    IntegerRaster susceptible;
    std::vector<cell> suitable_cells;
    std::vector<IntegerRaster> mortality_tracker_vector;
    IntegerRaster died;

    double establishment_probability_at(i, j)
    {
        double probability = double(this->susceptible_at(row, col))
                             / environment.total_population_at(row, col);
        if (weather)
            probability *= weather_coefficient(i, j);
        return probability;
    }
    double total_host_at(i, j)
    {
        return this->susceptible_at(row, col) + this->infected_at(row, col)
               + this->exposed_at(row, col);
    }
    double exposed_at(i, j)
    {
        sum = 0;
        for (const IntegerRaster& raster : exposed)
            sum += raster(i, j);
        return sum;
    }
    int susceptible_at(i, j)
    {
        return susceptible(i, j);
    }
    // interaction version
    // int susceptible_at(i, j)
    // {
    //     return environment.susceptible_at("Oak", i, j);
    // }
    void add_disperser_at(i, j)
    {
        if ("SI") {
            infected(i, j) += 1;
        }
        if ("SEI") {
            exposed.back()(i, j) += 1;
        }
    }
    void step_forward()  // aka infect
    {
        if ("SEI") {
            infected = exposed.front();
            rotate_left_by_one(exposed);
        }
        if (mortality) {
            rotate_left_by_one(mortality_tracker_vector);
        }
    }
};

template<typename IntegerRaster>
class Environment
{
    std::vector<Host*> hosts;  // non-owning reference
    IntegerRaster other_individuals;  // non hosts

    int total_population_at(i, j)
    {
        auto sum = other_individuals(i, j);
        for (const auto& host : hosts)
            sum += host.total_host_at(i, j);
        return sum;
    }
};

template<typename IntegerRaster>
class HostPool
{
    std::vector<Host> hosts;
    bool disperser_to(i, j)
    {
        std::vector probabilities;
        for (const auto& host : hosts) {
            probabilities.push_back(host.establishment_probability_at(i, j));
        }
        auto index = pick_host_index_by_probability(hosts, probabilities);
        probability = probabilities[index];
        if (establishment_stochasticity_)
            tester = distribution_uniform(generator_);
        else
            tester = 1 - deterministic_probability;
        if (tester < probability) {
            // establish
            host[index].add_disperser_at(i, j);
            return true;
        }
        // not established
        return false;
    }
};

template<typename IntegerRaster>
class PestPool
{
    IntegerRaster mobile_dispersers;
    IntegerRaster immobile_dispersers;  // aka in soil storage

    void immobile_dispersers_to(number, i, j)
    {
        // weather + ratios + stochasticity
        immobile_dispersers(i, j) += number;
    }
    int immobile_dispersers_from(i, j)
    {
        // weather + ratios + stochasticity
        number = 1;
        immobile_dispersers(i, j) -= number;
        return number;
    }
    int add_dispersers_at(dispersers, i, j)
    {
        // To soils
        dispersers_to_soil = to_soil_percentage * dispersers;
        this->immobile_dispersers_to(dispersers_to_soil, i, j);
        dispersers -= dispersers_to_soil;
        return dispersers;
    }
};

template<typename Hosts, typename Pests>
class Spread
{
    // Generate and kernel part of disperse could be the same,
    // then the disperser moves would be resolved right away,
    // and there would be no need to distinguish between generated dispersers,
    // but it would be much harder to track origin of established dispersers.

    void action(Hosts hosts, Pests pests)
    {
        for (const auto& [i, j] : hosts.suitable_cells) {
            // Generate
            dispersers = hosts.dispersers_from(i, j);
            // To soils
            dispersers = pests.add_dispersers_at(dispersers, i, j);

            // Classic spread
            for (const auto& disperser : dispersers) {
                i, j, kernel_name = kernel(disperser);
                // if (i, j is outside) {
                //    pests.add_free(i, j);
                //}
                pests.add_landed(i, j);
                established = hosts.disperser_to(i, j);
                if (established) {
                    pests.add_established(i, j, disperser, kernel_name);
                }
            }
            // From soils
            dispersers = pests.immobile_dispersers_from(i, j);
            for (const auto& disperser : dispersers) {
                established = hosts.disperser_to(i, j);
                if (established)
                    pests.add_established(i, j, i, j, "soil");
            }
        }
    }
};

template<typename Hosts>
class SurvivalRate
{
    void action(Hosts hosts)
    {
        for (const auto& [i, j] : hosts.suitable_cells) {
            count = hosts.infected_at(i, j) * survival_rate(i, j);
            hosts.remove_infection_at(i, j, count);
        }
    }
};

template<typename Hosts>
class RemoveByTemperature
{
    void action(Hosts hosts)
    {
        for (const auto& [i, j] : hosts.suitable_cells) {
            if (environment.temperature_at(i, j) < lethal_temperature) {
                count = hosts.infected_at(i, j);
                hosts.remove_infection_at(i, j, count);
            }
        }
    }
};

template<typename Hosts>
class PesticideTreatment
{
    void apply_treatment(Hosts hosts)
    {
        for (const auto& [i, j] : hosts.suitable_cells) {
            susceptible = treatment(i, j) * hosts.get_susceptible_at(i, j);
            infected = get_resistant(treatment(i, j), hosts.infected_at(i, j));
            vector<int> exposed_list;
            for
                value in hosts.exposed_at(i, j)
                {
                    exposed_list.push_back(get_resistant(treatment(i, j), value));
                }
            hosts.add_resistant_at(
                i,
                j,
                infected,
                susceptible,
                exposed_list);  // This could be an apply function with lambda.
        }
    }
    void void end_treatment_effect(Hosts hosts)
    {
        for (const auto& [i, j] : hosts.suitable_cells) {
            hosts.remove_all_resistant_at(i, j);
        }
    }
};

template<typename Hosts>
class Mortality
{
    void action(Hosts host)
    {
        // Do all mortality in Hosts and Host and here just call the right method.
        // Even mortality rate is linked to the host. Some don't die at all.
        for (const auto& [i, j] : hosts.suitable_cells) {
            for
                host in hosts
                {
                for
                    item in hosts.mortality_at(i, j)
                    {
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
};

template<typename Hosts>
class MoveOverpopulatedPests
{
    void action(Hosts hosts)
    {
        // generate_moves as from cell, to cell, count won't work because
        // 1) they need to be airborne first and
        // 2) only count for leaving=infected->susceptible is not tracked (no exposed
        // and resistant)
        for (const auto& [i, j] : hosts.suitable_cells) {
            infected = hosts.infected_at(i, j);
            if (infected / hosts.total_hosts_at(i, j) >= overpopulation_percentage) {
                i, j = kernel(disperser);
                hosts.disperser_from(i, j);
                moves.push_back(i, j);
            }
        }
        for
            i, j in moves
            {
                hosts.disperser_to(i, j);
            }
    }
};

/*
template<typename Hosts>
class Model
{
    Hosts hosts;
    hosts
        .set_infected(...)...  // more rasters added here (or all are set in the ctor)

        remove_by_temperature.action(hosts);
    survival_rate.action(hosts);

    move_overpopulated_pests.action(hosts);
    mortality.action(hosts);

    // The treatment would be wrapped in Treatments which has manage to decide
    // the time, but it does not need manage_mortality anymore because the
    // loop over mortality tracker is in Hosts.
    pesticide_treatment.apply_treatment(hosts);
    pesticide_treatment.end_treatment_effect(hosts);

    hosts.step_forward();
}
*/
}  // namespace pops
