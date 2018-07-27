#ifdef POPSS_TEST

/*
 * Simple compilation test for the PoPSS Raster class.
 *
 * Copyright (C) 2018 by the authors.
 *
 * Authors: Vaclav Petras <wenzeslaus gmail com>
 *
 * This file is part of PoPSS.

 * PoPSS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * PoPSS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with PoPSS. If not, see <https://www.gnu.org/licenses/>.
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


int main(int argc, char *argv[])
{
    Raster<int> a(10, 10, 1, 1);
    Raster<float> b(10, 10, 1, 1);
    Raster<double> c(10, 10, 1, 1);

    return 0;
}

#endif  // POPSS_TEST
