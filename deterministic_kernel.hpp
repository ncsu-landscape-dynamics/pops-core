/*
 * PoPS model - deterministic dispersal kernel
 *
 * Copyright (C) 2015-2020 by the authors.
 *
 * Authors: Margaret Lawrimore (malawrim ncsu edu)
 *
 * The code contained herein is licensed under the GNU General Public
 * License. You may obtain a copy of the GNU General Public License
 * Version 2 or later at the following locations:
 *
 * http://www.opensource.org/licenses/gpl-license.html
 * http://www.gnu.org/copyleft/gpl.html
 */

#ifndef POPS_DETERMINISTIC_KERNEL_HPP
#define POPS_DETERMINISTIC_KERNEL_HPP

#include "kernel_types.hpp"
#include "radial_kernel.hpp"

#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef PI
#define PI M_PI
#endif

namespace pops {

/*! 
 * Cauchy distribution
 * Includes probability density function and cumulative distribution function
 * pdf returns the probability that the variate has the value x
 * cdf returns the upper range that encompasses x percent of the distribution (e.g for 99% input .99)
 */
class cauchy_distribution
{
public:
	cauchy_distribution(double scale, double locator)
        : s(scale), t(locator)
    {}

    double pdf(double x)
    {
    	return 1 / ((s*PI)*(1 + (std::pow((x - t)/s, 2))));
    
    }
    
    double cdf(double x)
    {
    	if ( s == 1 && t == 0) {
    		// checking upper half so multiply result by 2
    		return std::tan((x / 2) * M_PI) * 2;
    	} else {
    		return (std::tan((x/2)*M_PI + (std::atan(-t/s))) * s + t) * 2;
    	}
    }
private:
	// scale parameter - 1 for standard
    double s;
    // location parameter - location of peak - 0 for standard
    double t;
};

/*! 
 * Exponential distribution
 * Includes probability density function and cumulative distribution function
 * pdf returns the probability that the variate has the value x
 * cdf returns the upper range that encompasses x percent of the distribution (e.g for 99% input 0.99)
 */
class exponential_distribution
{
public:
exponential_distribution(double scale)
        : beta(scale)
    {}
    // assumes mu is 0 which is traditionally accepted
    double pdf(double x)
    {
    	return (1/beta)*(std::exp(-x/beta));
    }
    
    double cdf(double x)
    {
    	if ( beta == 1) {
    		return -std::log(1-x);
    	} else {
    		return -beta * std::log(1-x);
    	}
    }
private:
	// scale parameter - 1 for standard 
	// equal to 1/lambda
	double beta;
};

/*! 
 * Dispersal kernel for deterministic spread to cell with highest probability of spread
 *
 * Dispersal Kernel type determines use of Exponential or Cauchy distribution
 * to find probability.
 *
 * dispersal_percentage is the percent of all possible dispersal to be included 
 * in the moving window size (e.g for 99% input 0.99).
 *
 * Useful for testing as it is deterministic and provides fully replicable results 
 */
class DeterministicDispersalKernel
{
protected:
    Raster<int> dispersers_;
    DispersalKernelType kernel_type_;
    // the west-east resolution of the pixel
    double east_west_resolution;
    // the north-south resolution of the pixel
    double north_south_resolution;
    // row/col position of middle cell
    int mid_row = 0;
    int mid_col = 0;
    // position of cell from previous call
    int prev_row = -1;
    int prev_col = -1;
    // number of rows/cols in the probability window
    int number_of_rows = 0;
    int number_of_columns = 0;
    // maximum distance from center cell to outer cells
    double max_distance;
    std::vector<std::vector<double>> probability;
    std::vector<std::vector<double>> probability_copy;
    cauchy_distribution cauchy;
    exponential_distribution exponential;
    double number_of_dispersers;
public:
    DeterministicDispersalKernel(DispersalKernelType dispersal_kernel, 
    							Raster<int> dispersers, double dispersal_percentage,
    							double ew_res, double ns_res, double distance_scale = 1,
    							double locator = 0)
        :
        dispersers_(dispersers),
        cauchy(distance_scale, locator),
        exponential(distance_scale),
        kernel_type_(dispersal_kernel),
        east_west_resolution(ew_res),
        north_south_resolution(ns_res)
    {
    	if (kernel_type_ == DispersalKernelType::Cauchy) {
    		max_distance = cauchy.cdf(dispersal_percentage);
    	} else if (kernel_type_ == DispersalKernelType::Exponential) {
    		max_distance = exponential.cdf(dispersal_percentage);
    	} else { 
    		//for expanding compatible kernel types
    	}
    	number_of_columns = std::ceil(max_distance / east_west_resolution) * 2 + 1;
    	number_of_rows = std::ceil(max_distance / north_south_resolution) * 2 + 1;
    	probability.resize(number_of_rows, std::vector<double>(number_of_columns));
    	probability_copy.resize(number_of_rows, std::vector<double>(number_of_columns));
    	mid_row = number_of_rows / 2;
    	mid_col = number_of_columns / 2;
    	double sum = 0;
    	for ( int i = 0; i < number_of_rows; i++ ) {
    		for ( int j = 0; j < number_of_columns; j++ ) {
    			double distance_to_center = 
    			std::sqrt(std::pow((std::abs(mid_row - i) * east_west_resolution), 2) 
    			+ pow((std::abs(mid_col - j) * north_south_resolution), 2));
    			// determine probability based on distance
				if (kernel_type_ == DispersalKernelType::Cauchy) {
 					probability[i][j] = std::abs(cauchy.pdf(distance_to_center));
 				} else if ( kernel_type_ == DispersalKernelType::Exponential ) {
 					probability[i][j] = std::abs(exponential.pdf(distance_to_center));
				}
 				sum += probability[i][j];
    		}
    	}
    	//normalize based on the sum of all probabilities in the raster
    	for ( int i = 0; i < number_of_rows; i++ ) {
    		for ( int j = 0; j < number_of_columns; j++ ) {
    			probability[i][j] /= sum;
    		}
    	}
    }

    /*! Generates a new position for the spread.
     *
     *  Creates a copy of the probability matrix to mark where dispersers are assigned.
     *  New window created any time a new cell is selected from simulation.disperse
     *
     *  Selects next row/col value based on the cell with the highest probability
     *  in the window.
     *
     *  Generator is not used
     */
    template<typename Generator>
    std::tuple<int, int> operator() (Generator& generator, int row, int col)
    {
		// reset the window if considering a new cell
    	if ( row != prev_row || col != prev_col ) {
    		number_of_dispersers = 1.0 / dispersers_(row, col);
    		for ( int i = 0; i < number_of_rows; i++) {
    			for ( int j = 0; j < number_of_columns; j++) {
    				probability_copy[i][j] = probability[i][j];
    			}
    		}
    	}
    	
    	int row_movement = 0;
    	int col_movement = 0;
    	
    	double max = (double)-std::numeric_limits<int>::max();
    	int max_prob_row = 0;
    	int max_prob_col = 0;
    	
    	//find cell with highest probability
    	for ( int i = 0; i < number_of_rows; i++ ) {
    		for ( int j = 0; j < number_of_columns; j++ ) {
    			if (probability_copy[i][j] > max ) {
    				max = probability_copy[i][j];
    				max_prob_row = i;
					max_prob_col = j;
					row_movement = i - mid_row;
					col_movement = j - mid_col;
    			}
    		}
    	}

		// subtracting 1/number of dispersers ensures we always move the same proportion 
		// of the individuals to each cell no matter how many are dispersing
    	probability_copy[max_prob_row][max_prob_col] -= number_of_dispersers;
    	prev_row = row;
    	prev_col = col;
    	
    	// need to return values in terms of actual location
        return std::make_tuple(row + row_movement, col + col_movement);
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
    static bool supports_kernel(const DispersalKernelType type)
    {
        static const std::array<DispersalKernelType, 2> supports = {
            DispersalKernelType::Cauchy,
            DispersalKernelType::Exponential
        };
        auto it = std::find(supports.cbegin(), supports.cend(), type);
        return it != supports.cend();;
    }
};

} // namespace pops

#endif // POPS_DETERMINISTIC_KERNEL_HPP
