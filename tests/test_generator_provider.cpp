/*
 * Test random number generator provider classes.
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

#include <random>

#include <pops/generator_provider.hpp>

using namespace pops;

/**
 * Return true if the *operation* throws an exception.
 *
 * Specify exception using the template parameter and the operation using a lambda
 * function.
 *
 */
template<typename Exception, typename Operation>
bool throws_exception(Operation operation)
{
    try {
        operation();
    }
    catch (const Exception&) {
        return true;
    }
    return false;
}

/** Return 1 and print message if two number are different, zero otherwise */
int assert_pair_equals(
    std::string test_name,
    int test_number,
    int first,
    int second,
    std::string first_name,
    std::string second_name)
{
    if (first != second) {
        std::cerr << test_name << " (" << test_number << "): " << first_name << ": "
                  << first << ", " << second_name << ": " << second
                  << " should equal, but are different\n";
        return 1;
    }
    return 0;
}

/** Return 1 and print message if two number are the same, zero otherwise */
int assert_pair_not_equals(
    std::string test_name,
    int test_number,
    int first,
    int second,
    std::string first_name,
    std::string second_name)
{
    if (first == second) {
        std::cerr << test_name << " (" << test_number << "): " << first_name << ": "
                  << first << ", " << second_name << ": " << second
                  << " should be different, but are equal\n";
        return 1;
    }
    return 0;
}

/** Return 1 and print message if values don't equal, 0 otherwise. */
template<typename Type>
int assert_value_equals_value(std::string test_name, Type first, Type second)
{
    if (first != second) {
        std::cerr << test_name << ": " << first << " != " << second << "\n";
        return 1;
    }
    return 0;
}

/** Return 1 and print message if key is not present or value is not expected, 0
 * otherwise. */
int assert_value_under_key(
    std::string test_name,
    const std::map<std::string, unsigned>& seeds,
    std::string key,
    unsigned value)
{
    int ret = 0;
    if (!container_contains(seeds, key)) {
        std::cerr << test_name << ": key '" << key << "' not in the map\n";
        ++ret;
    }
    if (!ret)
        ret += assert_value_equals_value<unsigned>(test_name, seeds.at(key), value);
    if (ret) {
        std::cerr << "map contains:\n";
        for (const auto& item : seeds) {
            std::cerr << " " << item.first << ": " << item.second << "\n";
        }
    }
    return ret;
}

/**
 * Check that accessing different simple generators gives the same result as accesing
 * a single one when the seed is the same.
 */
int test_single_generator_results_same()
{
    int ret = 0;
    unsigned seed = 42;
    DefaultSingleGeneratorProvider generator1(seed);
    DefaultSingleGeneratorProvider generator2(seed);
    int a = 13;
    int b = 27;
    std::uniform_int_distribution<int> distribution1(a, b);
    std::uniform_int_distribution<int> distribution2(a, b);
    int repetions = 10;
    for (int i = 0; i < repetions; ++i) {
        int number1 = distribution1(generator1.weather());
        int number2 = distribution2(generator2.disperser_generation());
        ret += assert_pair_equals(
            "test_single_generator_results_same",
            i,
            number1,
            number2,
            "weather",
            "disperser_generation");
        number1 = distribution1(generator1.establishment());
        number2 = distribution2(generator2.disperser_generation());
        ret += assert_pair_equals(
            "test_single_generator_results_same",
            i,
            number1,
            number2,
            "establishment()",
            "disperser_generation");
        number1 = distribution1(generator1.movement());
        number2 = distribution2(generator2.disperser_generation());
        ret += assert_pair_equals(
            "test_single_generator_results_same",
            i,
            number1,
            number2,
            "movement",
            "disperser_generation");
    }
    return ret;
}

/**
 * Check that accessing different generators gives the same result as accesing
 * a single one when the initial seed is the same but independent seeds are used
 * within each provider.
 */
int test_multiple_generator_results_same()
{
    int ret = 0;
    unsigned seed = 42;
    RandomNumberGeneratorProvider<std::default_random_engine> generator1(seed, true);
    RandomNumberGeneratorProvider<std::default_random_engine> generator2(seed, true);
    int a = 13;
    int b = 27;
    std::uniform_int_distribution<int> distribution1(a, b);
    std::uniform_int_distribution<int> distribution2(a, b);
    int repetions = 10;
    for (int i = 0; i < repetions; ++i) {
        int number1 = distribution1(generator1.weather());
        int number2 = distribution2(generator2.weather());
        ret += assert_pair_equals(
            "test_multiple_generator_results_same",
            i,
            number1,
            number2,
            "weather 1",
            "weather 2 ");
        number1 = distribution1(generator1.overpopulation());
        number2 = distribution1(generator2.overpopulation());
        ret += assert_pair_equals(
            "test_multiple_generator_results_same",
            i,
            number1,
            number2,
            "overpopulation 1",
            "overpopulation 2");
    }
    return ret;
}

/**
 * Check that accessing different generators gives the same result as accesing
 * a single one when the seed is the same when one is used differently.
 * (This is similar to how different model runs would turn off and on different parts
 * of the model.)
 */
int test_multiple_generator_results_independent()
{
    int ret = 0;
    unsigned seed = 42;
    RandomNumberGeneratorProvider<std::default_random_engine> generator1(seed, true);
    RandomNumberGeneratorProvider<std::default_random_engine> generator2(seed, true);
    int a = 13;
    int b = 27;
    std::uniform_int_distribution<int> distribution1(a, b);
    std::uniform_int_distribution<int> distribution2(a, b);
    int repetions = 10;
    for (int i = 0; i < repetions; ++i) {
        int number1 = distribution1(generator1.weather());
        int number2 = distribution2(generator2.weather());
        ret += assert_pair_equals(
            "test_multiple_generator_results_independent - first attempt",
            i,
            number1,
            number2,
            "generator 1",
            "generator 2");
        // Use another generator, but only from one provider.
        int number3 = distribution1(generator1.overpopulation());
        number1 = distribution1(generator1.weather());
        number2 = distribution2(generator2.weather());
        ret += assert_pair_equals(
            "test_multiple_generator_results_independent "
            "- second attempt after drawing "
                + std::to_string(number3) + " from overpopulation",
            i,
            number1,
            number2,
            "generator 1",
            "generator 2");
    }
    return ret;
}

/**
 * Check that multiple custom seeds give same results for same seeds and different
 * results for different seeds.
 */
int test_multiple_seeds()
{
    int ret = 0;
    std::map<std::string, unsigned> seeds{
        {{"disperser_generation", 42},
         {"natural_dispersal", 342},
         {"anthropogenic_dispersal", 342},
         {"establishment", 462},
         {"weather", 252},
         {"lethal_temperature", 101},
         {"movement", 72},
         {"overpopulation", 42},
         {"survival_rate", 252},
         {"soil", 462}}};
    RandomNumberGeneratorProvider<std::default_random_engine> generator(seeds);
    int a = 13;
    int b = 1278;  // Wide range to minimize overlap by chance for some seeds.
    std::uniform_int_distribution<int> disperser_generation_distribution(a, b);
    std::uniform_int_distribution<int> weather_distribution(a, b);
    std::uniform_int_distribution<int> establishment_distribution(a, b);
    std::uniform_int_distribution<int> movement_distribution(a, b);
    std::uniform_int_distribution<int> overpopulation_distribution(a, b);
    std::uniform_int_distribution<int> survival_rate_distribution(a, b);
    std::uniform_int_distribution<int> soil_distribution(a, b);
    int repetions = 10;
    for (int i = 0; i < repetions; ++i) {
        int disperser_generation =
            disperser_generation_distribution(generator.disperser_generation());
        int weather = weather_distribution(generator.weather());
        int establishment = establishment_distribution(generator.establishment());
        int movement = movement_distribution(generator.movement());
        int overpopulation = overpopulation_distribution(generator.overpopulation());
        int survival_rate = survival_rate_distribution(generator.survival_rate());
        int soil = soil_distribution(generator.soil());
        ret += assert_pair_equals(
            "test_multiple_seeds",
            i,
            disperser_generation,
            overpopulation,
            "disperser_generation",
            "overpopulation");
        ret += assert_pair_equals(
            "test_multiple_seeds",
            i,
            weather,
            survival_rate,
            "weather",
            "survival_rate");
        ret += assert_pair_equals(
            "test_multiple_seeds", i, establishment, soil, "establishment", "soil");
        // These are seed-dependent.
        ret += assert_pair_not_equals(
            "test_multiple_seeds",
            i,
            establishment,
            movement,
            "establishment",
            "movement");
        ret += assert_pair_not_equals(
            "test_multiple_seeds",
            i,
            overpopulation,
            survival_rate,
            "overpopulation",
            "survival_rate");
        ret += assert_pair_not_equals(
            "test_multiple_seeds", i, survival_rate, soil, "survival_rate", "soil");
    }
    return ret;
}

/** Check that reading of seeds works for key=value parameter style */
int test_seed_config_parameter_style()
{
    int ret = 0;
    std::string text("weather = 252 , establishment =562,survival_rate=252,soil=462");
    Config config;
    config.read_seeds(text, ',', '=');
    const auto& seeds = config.random_seeds;
    ret += assert_value_under_key(
        "test_seed_config_parameter_style", seeds, "weather", 252);
    ret += assert_value_under_key(
        "test_seed_config_parameter_style", seeds, "establishment", 562);
    ret += assert_value_under_key(
        "test_seed_config_parameter_style", seeds, "survival_rate", 252);
    ret +=
        assert_value_under_key("test_seed_config_parameter_style", seeds, "soil", 462);
    if (config.multiple_random_seeds != true) {
        std::cerr << "test_seed_config_parameter_style: config.multiple_random_seeds"
                     " not set to true\n";
        ++ret;
    }
    return ret;
}

/** Check that reading of seeds works for YAML config style */
int test_seed_config_yaml_style()
{
    int ret = 0;
    std::string text("weather: 252\nestablishment: 562\nsurvival_rate:252\nsoil:462");
    Config config;
    config.read_seeds(text, '\n', ':');
    const auto& seeds = config.random_seeds;
    ret += assert_value_under_key("test_seed_config_yaml_style", seeds, "weather", 252);
    ret += assert_value_under_key(
        "test_seed_config_yaml_style", seeds, "establishment", 562);
    ret += assert_value_under_key(
        "test_seed_config_yaml_style", seeds, "survival_rate", 252);
    ret += assert_value_under_key("test_seed_config_yaml_style", seeds, "soil", 462);
    if (config.multiple_random_seeds != true) {
        std::cerr << "test_seed_config_parameter_style: config.multiple_random_seeds"
                     " not set to true\n";
        ++ret;
    }
    return ret;
}

/**
 * Check that ordered seeds are accepted and all are assigned
 *
 * Importantly, this tests that the two parts of source code which have
 * a list of seeds by name are in sync.
 */
int test_seed_config_vector_works()
{
    int ret = 0;
    Config config;
    // Letting exceptions handle the error states here.
    config.read_seeds({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    bool thrown = throws_exception<std::invalid_argument>(
        [&config] { validate_random_number_generator_provider_config(config); });
    if (thrown) {
        std::cerr << "test_seed_config_vector_works: "
                     "List of seed names in Config is not in sync with providers\n";
        ++ret;
    }
    return ret;
}

/** Check that reading of seeds from list fails if there is not enough values */
int test_seed_config_vector_incomplete_list_recognized()
{
    int ret = 0;
    Config config;
    bool thrown = throws_exception<std::invalid_argument>(
        [&config] { config.read_seeds({1, 2, 3}); });
    if (!thrown) {
        std::cerr << "test_seed_config_vector_incomplete_list_recognized: "
                     "An incomplete list of seeds wrongly accepted\n";
        ++ret;
    }
    validate_random_number_generator_provider_config(config);
    return ret;
}

/** Check that a valid list of seeds successfully validates */
int test_seed_config_validate_validates_valid()
{
    std::string text(
        "disperser_generation=1,natural_dispersal=2,anthropogenic_dispersal=1,"
        "establishment=2,weather=2,lethal_temperature=20,establishment=3,movement=4,"
        "overpopulation=5,survival_rate=6,soil=7");
    Config config;
    config.read_seeds(text, ',', '=');
    validate_random_number_generator_provider_config(config);
    validate_random_number_generator_provider_seeds(config.random_seeds);
    // Letting exceptions handle the error states.
    return 0;
}

/** Check that incomplete list of seeds causes exception */
int test_seed_config_validate_throws()
{
    int ret = 0;
    std::string text("weather=252,establishment=562,survival_rate=252,soil=462");
    Config config;
    config.read_seeds(text, ',', '=');
    bool thrown = throws_exception<std::invalid_argument>(
        [&config] { validate_random_number_generator_provider_config(config); });
    if (!thrown) {
        std::cerr << "An incomplete list of seeds wrongly accepted "
                  << "by validate_random_number_generator_provider_config()\n";
        ++ret;
    }
    thrown = throws_exception<std::invalid_argument>([&config] {
        validate_random_number_generator_provider_seeds(config.random_seeds);
    });
    if (!thrown) {
        std::cerr << "An incomplete list of seeds wrongly accepted "
                  << "by validate_random_number_generator_provider_seeds()\n";
        ++ret;
    }
    return ret;
}

/** Run all tests and collect the resulting return value. */
int run_tests()
{
    int ret = 0;

    ret += test_single_generator_results_same();
    ret += test_multiple_generator_results_same();
    ret += test_multiple_generator_results_independent();
    ret += test_multiple_seeds();
    ret += test_seed_config_parameter_style();
    ret += test_seed_config_yaml_style();
    ret += test_seed_config_vector_works();
    ret += test_seed_config_vector_incomplete_list_recognized();
    ret += test_seed_config_validate_validates_valid();
    ret += test_seed_config_validate_throws();

    if (ret)
        std::cerr << "Number of errors in the generator provider test: " << ret << "\n";
    return ret;
}

/** Run the test or error when command line parameters are provided. */
int main(int argc, char**)
{
    if (argc > 1)
        return 1;
    else
        return run_tests();
}
