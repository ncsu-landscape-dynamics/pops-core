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
    SingleGeneratorProvider generator1(seed);
    SingleGeneratorProvider generator2(seed);
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

int run_tests()
{
    int ret = 0;

    ret += test_generator();

    if (ret)
        std::cerr << "Number of errors in the network test: " << ret << "\n";
    return ret;
}

int main(int argc, char** argv)
{
    if (argc > 1)
        return 1;
    else
        return run_tests();
}
