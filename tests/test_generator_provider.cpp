#include <fstream>
#include <regex>
#include <random>

#include <pops/simple_generator.hpp>

using namespace pops;

/**
 * Check that accessing different generators gives the same result as accesing
 * a single one.
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
        int number2 = distribution2(generator2.general());
        if (number1 != number2) {
            std::cerr << "test simple generator - weather (" << i
                      << "): generator 1: " << number1 << " generator 2: " << number2
                      << "\n";
            ret += 1;
        }
        number1 = distribution1(generator1.lethal_temperature());
        number2 = distribution2(generator2.general());
        if (number1 != number2) {
            std::cerr << "test simple generator - lethal_temperature (" << i
                      << "): generator 1: " << number1 << " generator 2: " << number2
                      << "\n";
            ret += 1;
        }
        number1 = distribution1(generator1.movement());
        number2 = distribution2(generator2.general());
        if (number1 != number2) {
            std::cerr << "test simple generator - movement() (" << i
                      << "): generator 1: " << number1 << " generator 2: " << number2
                      << "\n";
            ret += 1;
        }
    }
    return ret;
}

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
        if (number1 != number2) {
            std::cerr << "test multiple generator - weather (" << i
                      << "): generator 1: " << number1 << " generator 2: " << number2
                      << "\n";
            ret += 1;
        }
        number1 = distribution1(generator1.overpopulation());
        number2 = distribution1(generator2.overpopulation());
        if (number1 != number2) {
            std::cerr << "test multiple generator - overpopulation (" << i
                      << "): generator 1: " << number1 << " generator 2: " << number2
                      << "\n";
            ret += 1;
        }
    }
    return ret;
}

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
        if (number1 != number2) {
            std::cerr << "test multiple generator - weather (" << i
                      << "): generator 1: " << number1 << " generator 2: " << number2
                      << "\n";
            ret += 1;
        }
        // Use another generator, but only from one provider.
        int number3 = distribution1(generator1.overpopulation());
        number1 = distribution1(generator1.weather());
        number2 = distribution2(generator2.weather());
        if (number1 != number2) {
            std::cerr << "test multiple generator - weather - second attempt (" << i
                      << "): generator 1: " << number1 << " generator 2: " << number2
                      << "(another number was: " << number3 << "\n";
            ret += 1;
        }
    }
    return ret;
}

int run_tests()
{
    int ret = 0;

    ret += test_single_generator_results_same();
    ret += test_multiple_generator_results_independent();
    ret += test_multiple_generator_results_independent();

    if (ret)
        std::cerr << "Number of errors in the network test: " << ret << "\n";
    return ret;
}

int main(int argc, char**)
{
    if (argc > 1)
        return 1;
    else
        return run_tests();
}
