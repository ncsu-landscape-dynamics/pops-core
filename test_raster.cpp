#ifdef POPS_TEST

/*
 * Simple compilation test for the PoPS Raster class.
 *
 * Copyright (C) 2018 by the authors.
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

#include "raster.hpp"
#include "simulation.hpp"

#include <map>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>

using std::string;
using std::cout;
using std::cerr;
using std::endl;

using pops::Raster;

static
void test_constructor_by_type()
{
    Raster<int> a(10, 10);
    Raster<float> b(10, 10);
    Raster<double> c(10, 10);
}

static
void test_constructor_dimensions()
{
    int x = 5;
    int y = 3;
    // new based on dimensions, with zeros, one value in the corner
    Raster<int> a(x, y);
    a.zero();
    a(x - 1, y - 1) = 2;
    std::cout << x << "x" << y << ":" << std::endl << a;
    // same as above, but with the dimensions swapped
    Raster<int> b(y, x);
    b.zero();
    b(y - 1, x - 1) = 2;
    std::cout << y << "x" << x << ":" << std::endl << b;
}

static
void test_initializer_and_output()
{
    Raster<int> a = {{1, 2},
                     {3, 4},
                     {5, 6}};
    std::cout << "3x2:" << std::endl << a;

    Raster<int> b = {{0, 0, 0, 5, 5},
                     {0, 0, 0, 5, 5},
                     {0, 0, 0, 10, 0}};
    std::cout << "3x5:" << std::endl << b;
}

static
void test_equal_operator()
{
    Raster<int> a = {{1, 2}, {3, 4}, {5, 6}};
    Raster<int> b = {{1, 2}, {3, 4}, {5, 6}};
    if (a == b)
        std::cout << "Operator equal works" << std::endl;
    else
        std::cout << "Operator equal does not work" << std::endl;
}

static
void test_not_equal_operator()
{
    Raster<int> a = {{1, 2}, {3, 4}, {5, 6}};
    Raster<int> b = {{1, 2}, {3, 5}, {5, 6}};
    if (a != b)
        std::cout << "Operator not-equal works" << std::endl;
    else
        std::cout << "Operator not-equal does not work" << std::endl;
}

static
int test_plus_operator()
{
    Raster<int> d = {{1, 2}, {3, 4}, {5, 6}};
    Raster<int> e = {{8, 9}, {10, 11}, {12, 13}};
    Raster<int> f = {{9, 11}, {13, 15}, {17, 19}};
    if (d + e == f) {
        std::cout << "Operator plus works" << std::endl;
        return 0;
    }
    else {
        std::cout << "Operator plus does not work" << std::endl;
        return 1;
    }
}

static
int test_multiply_in_place_operator()
{
    Raster<double> d = {{1.1, 2}, {3.84, 4}, {5, 6}};
    Raster<double> e = {{8, 9.5}, {10, 11}, {12, 13}};
    Raster<double> f = e;
    Raster<double> g = {{8.8, 19}, {38.4, 44}, {60, 78}};
    d *= e;
    if (e == f && d == g) {
        std::cout << "Operator *= works" << std::endl;
        return 0;
    }
    else {
        std::cout << "Operator *= does not work" << std::endl;
        return 1;
    }
}

static
void test_sqrt()
{
    Raster<int> a = {{16, 25}, {4, 9}};
    Raster<int> b = {{4, 5}, {2, 3}};
    auto c = sqrt(a);
    std::cout << "sqrt function: ";
    if (b == c)
        std::cout << "OK" << std::endl;
    else
        std::cout << "\n" << a << "!=\n" << b << std::endl;
}

template<typename T, typename U, typename V>
static
void test_diff_types_const()
{
    Raster<T> a = {{1, 2}, {3, 4}, {5, 6}};
    Raster<U> b = {{1, 2}, {3, 4}, {5, 6}};
    Raster<V> c = a + b;
    c = a - b;
    c = a * b;
    c = a / b;
}

template<typename T, typename U>
static
void test_diff_types_modify()
{
    Raster<T> a = {{1, 2}, {3, 4}, {5, 6}};
    Raster<U> b = {{1, 2}, {3, 4}, {5, 6}};
    a += b;
    a -= b;
    a *= b;
    a /= b;
}

template<typename T>
static
void test_op_int()
{
    Raster<T> a = {{1, 2}, {3, 4}, {5, 6}};
    a += 1;
    a -= 10;
    a *= 6;
    a /= 2;
}

template<typename T>
static
void test_op_double()
{
    Raster<T> a = {{1, 2}, {3, 4}, {5, 6}};
    a += 1.1;
    a -= 10.3;
    a *= 7.5;
    a /= 2.5;
}

template<typename T>
static
void test_op_order()
{
    Raster<T> a = {{1, 2}, {3, 4}, {5, 6}};
    a + 1.1;
    a - 10.3;
    a * 75;
    a / 25;
    1 + a;
    10 - a;
    6.2 * a;
    2.1 / a;
}

int main()
{
    test_constructor_by_type();
    test_constructor_dimensions();
    test_initializer_and_output();
    test_equal_operator();
    test_not_equal_operator();

    test_plus_operator();
    test_multiply_in_place_operator();

    test_sqrt();

    // all doubles, no problem
    test_diff_types_const<double, double, double>();
    // operation on ints gives ints
    test_diff_types_const<int, int, int>();
    // combining int and double gives double
    test_diff_types_const<double, int, double>();
    test_diff_types_const<int, double, double>();
    // undefined: we don't define promotion
    // test_diff_types_const<int, int, double>();
    // undefined: we don't allow demotion
    // test_diff_types_const<double, double, int>();

    // all doubles, no problem
    test_diff_types_modify<double, double>();
    // operation on ints gives ints
    test_diff_types_modify<int, int>();
    // operation on double gives double
    test_diff_types_modify<double, int>();
    // undefined: modifying int with double not allowed
    // test_diff_types_modify<int, double>();

    test_op_double<double>();
    test_op_int<double>();
    test_op_int<int>();
    test_op_double<int>();

    test_op_order<double>();
    test_op_order<int>();

    return 0;
}

#endif  // POPS_TEST
