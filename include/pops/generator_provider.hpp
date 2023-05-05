/*
 * Simple random number generator provider class.
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

#ifndef POPS_SIMPLE_GENERATOR_HPP
#define POPS_SIMPLE_GENERATOR_HPP

#include <random>

#include "config.hpp"

namespace pops {

template<typename Generator>
class RandomNumberGeneratorProviderInterface
{
public:
    virtual void seed(unsigned seed) = 0;
    virtual void seed(const std::map<std::string, unsigned>& seeds) = 0;
    virtual void seed(const Config& config) = 0;
    virtual Generator& disperser_generation() = 0;
    virtual Generator& natural_dispersal() = 0;
    virtual Generator& anthropogenic_dispersal() = 0;
    virtual Generator& establishment() = 0;
    virtual Generator& weather() = 0;
    virtual Generator& lethal_temperature() = 0;
    virtual Generator& movement() = 0;
    virtual Generator& overpopulation() = 0;
    virtual Generator& survival_rate() = 0;
    virtual Generator& soil() = 0;
    virtual ~RandomNumberGeneratorProviderInterface() = default;
};

template<typename Generator>
class SingleGeneratorProvider : public RandomNumberGeneratorProviderInterface<Generator>
{
public:
    /**
     * @brief RandomNumberGeneratorProvider
     * @param seed
     *
     * Seeds first generator with the seed and then each subsequent generator with
     * seed += 1.
     */
    SingleGeneratorProvider(unsigned seed)
    {
        this->seed(seed);
    }

    void seed(unsigned seed)
    {
        general_generator_.seed(seed);
    }

    void seed(const std::map<std::string, unsigned>& seeds)
    {
        UNUSED(seeds);
        throw std::invalid_argument(
            "Multiple seeds are not supported by SimpleGeneratorProvider (only one seed is supported)");
    }

    void seed(const Config& config)
    {
        if (config.multiple_random_seeds) {
            throw std::invalid_argument(
                "Config cannot have multiple_random_seeds set for SimpleGeneratorProvider (only random_seed is supported)");
        }
        if (!config.random_seeds.empty()) {
            throw std::invalid_argument(
                "Config cannot have random_seeds set for SimpleGeneratorProvider (only random_seed is supported)");
        }
        general_generator_.seed(config.random_seed);
    }

    Generator& general()
    {
        return general_generator_;
    }

    Generator& disperser_generation()
    {
        return general();
    }

    Generator& natural_dispersal()
    {
        return general();
    }

    Generator& anthropogenic_dispersal()
    {
        return general();
    }

    Generator& establishment()
    {
        return general();
    }

    Generator& weather()
    {
        return general();
    }

    Generator& lethal_temperature()
    {
        return general();
    }

    Generator& movement()
    {
        return general();
    }

    Generator& overpopulation()
    {
        return general();
    }

    Generator& survival_rate()
    {
        return general();
    }

    Generator& soil()
    {
        return general();
    }

    // API to behave like the underlying generator.

    using result_type = typename Generator::result_type;

    static result_type min()
    {
        return Generator::min();
    }
    static result_type max()
    {
        return Generator::max();
    }
    result_type operator()()
    {
        return general_generator_();
    }

    void discard(unsigned long long n)
    {
        general_generator_.disard(n);
    }

private:
    Generator general_generator_;
};

using DefaultSingleGeneratorProvider =
    SingleGeneratorProvider<std::default_random_engine>;

template<typename Generator>
class IsolatedRandomNumberGeneratorProvider
    : public RandomNumberGeneratorProviderInterface<Generator>
{
public:
    /**
     * @brief RandomNumberGeneratorProvider
     * @param seed
     *
     * Seeds first generator with the seed and then each subsequent generator with
     * seed += 1.
     */
    IsolatedRandomNumberGeneratorProvider(unsigned seed)
    {
        this->seed(seed);
    }
    IsolatedRandomNumberGeneratorProvider(const std::map<std::string, unsigned>& seeds)
    {
        this->seed(seeds);
    }

    IsolatedRandomNumberGeneratorProvider(const Config& config)
    {
        this->seed(config);
    }

    void seed(unsigned seed)
    {
        disperser_generation_generator_.seed(seed++);
        natural_dispersal_generator_.seed(seed++);
        anthropogenic_dispersal_generator_.seed(seed++);
        establishment_generator_.seed(seed++);
        weather_generator_.seed(seed++);
        lethal_temperature_generator_.seed(seed++);
        movement_generator_.seed(seed++);
        overpopulation_generator_.seed(seed++);
        survival_rate_generator_.seed(seed++);
        soil_generator_.seed(seed);
    }

    void seed(const std::map<std::string, unsigned>& seeds)
    {
        this->set_seed_by_name(
            seeds, "disperser_generation", disperser_generation_generator_);
        this->set_seed_by_name(
            seeds, "natural_dispersal", natural_dispersal_generator_);
        this->set_seed_by_name(
            seeds, "anthropogenic_dispersal", anthropogenic_dispersal_generator_);
        this->set_seed_by_name(seeds, "establishment", establishment_generator_);
        this->set_seed_by_name(seeds, "weather", weather_generator_);
        this->set_seed_by_name(
            seeds, "lethal_temperature", lethal_temperature_generator_);
        this->set_seed_by_name(seeds, "movement", movement_generator_);
        this->set_seed_by_name(seeds, "overpopulation", overpopulation_generator_);
        this->set_seed_by_name(seeds, "survival_rate", survival_rate_generator_);
        this->set_seed_by_name(seeds, "soil", soil_generator_);
    }

    void seed(const Config& config)
    {
        if (!config.random_seeds.empty()) {
            this->seed(config.random_seeds);
        }
        else {
            this->seed(config.random_seed);
        }
    }

    Generator& disperser_generation()
    {
        return disperser_generation_generator_;
    }

    Generator& natural_dispersal()
    {
        return natural_dispersal_generator_;
    }

    Generator& anthropogenic_dispersal()
    {
        return anthropogenic_dispersal_generator_;
    }

    Generator& establishment()
    {
        return establishment_generator_;
    }

    Generator& weather()
    {
        return weather_generator_;
    }

    Generator& lethal_temperature()
    {
        return lethal_temperature_generator_;
    }

    Generator& movement()
    {
        return movement_generator_;
    }

    Generator& overpopulation()
    {
        return overpopulation_generator_;
    }

    Generator& survival_rate()
    {
        return survival_rate_generator_;
    }

    Generator& soil()
    {
        return soil_generator_;
    }

private:
    void set_seed_by_name(
        const std::map<std::string, unsigned>& seeds,
        const char* key,
        Generator& generator)
    {
        try {
            generator.seed(seeds.at(key));
        }
        catch (const std::out_of_range&) {
            throw std::invalid_argument(
                std::string("Seed '") + key + "' missing from the seeds configuration");
        }
    }

    Generator disperser_generation_generator_;
    Generator natural_dispersal_generator_;
    Generator anthropogenic_dispersal_generator_;
    Generator establishment_generator_;
    Generator weather_generator_;
    Generator lethal_temperature_generator_;  // Not need at this point.
    Generator movement_generator_;
    Generator overpopulation_generator_;
    Generator survival_rate_generator_;
    Generator soil_generator_;
};

template<typename Generator>
class RandomNumberGeneratorProvider
{
public:
    /**
     * @brief RandomNumberGeneratorProvider
     * @param seed
     *
     * Seeds first generator with the seed and then each subsequent generator with
     * seed += 1.
     */
    RandomNumberGeneratorProvider(unsigned seed, bool isolated = false)
        : isolated_(isolated)
    {
        if (isolated) {
            impl.reset(new IsolatedRandomNumberGeneratorProvider<Generator>(seed));
        }
        else {
            impl.reset(new SingleGeneratorProvider<Generator>(seed));
        }
    }
    RandomNumberGeneratorProvider(const std::map<std::string, unsigned>& seeds)
        : impl(new IsolatedRandomNumberGeneratorProvider<Generator>(seeds)),
          isolated_(true)
    {}

    RandomNumberGeneratorProvider(const Config& config) : impl(nullptr)
    {
        if (config.multiple_random_seeds) {
            impl.reset(new IsolatedRandomNumberGeneratorProvider<Generator>(config));
            isolated_ = true;
        }
        else {
            impl.reset(new SingleGeneratorProvider<Generator>(config.random_seed));
            isolated_ = false;
        }
    }

    void seed(unsigned seed)
    {
        return impl->seed(seed);
    }

    Generator& disperser_generation()
    {
        return impl->disperser_generation();
    }

    Generator& natural_dispersal()
    {
        return impl->natural_dispersal();
    }

    Generator& anthropogenic_dispersal()
    {
        return impl->anthropogenic_dispersal();
    }

    Generator& establishment()
    {
        return impl->establishment();
    }

    Generator& weather()
    {
        return impl->weather();
    }

    Generator& lethal_temperature()
    {
        return impl->lethal_temperature();
    }

    Generator& movement()
    {
        return impl->movement();
    }

    Generator& overpopulation()
    {
        return impl->overpopulation();
    }

    Generator& survival_rate()
    {
        return impl->survival_rate();
    }

    Generator& soil()
    {
        return impl->soil();
    }

    // API to behave like the underlying generator.

    using result_type = typename Generator::result_type;

    static result_type min()
    {
        return Generator::min();
    }
    static result_type max()
    {
        return Generator::max();
    }
    result_type operator()()
    {
        if (isolated_) {
            std::runtime_error(
                "RandomNumberGeneratorProvider used as a single generator "
                "but it is set to provide isolated generators");
        }
        return impl->disperser_generation().operator()();
    }

    void discard(unsigned long long n)
    {
        if (isolated_) {
            std::runtime_error(
                "RandomNumberGeneratorProvider used as a single generator "
                "but it is set to provide isolated generators");
        }
        impl->disperser_generation().disard(n);
    }

private:
    std::unique_ptr<RandomNumberGeneratorProviderInterface<Generator>> impl;
    bool isolated_ = false;
};

}  // namespace pops

#endif  // POPS_SIMPLE_GENERATOR_HPP
