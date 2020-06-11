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
 * includes probability density function and cumulative distribution function
 * pdf returns the probability that the variate has the value x
 * cdf returns the range that encompasses x percent of the distribution
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
    
    // x is the percentage of the distribution you want covered (e.g for 99% input .99)
    // rewrote cdf to not assume the integral is equal to 1, instead made it equal to 
    // x and solved for the maximum value which should return the value on the x-axis that
    // corresponds to the upper (x/2 %) of the distribution
    double cdf(double x)
    {
    	if ( s == 1 && t == 0) {
    		return std::tan((x / 2) * M_PI) * 2;
    		// need to multiply result by 2 since checking just the upper half of the
    		// symmetrical distribution
    	} else {
    		return (std::tan((x/2)*M_PI + (std::atan(-t/s))) * s + t) * 2;
    	}
    }
private:
	// scale parameter - specifies half width at half maximum - 1 for standard
	// often referred to as lamba which equals 1/beta
    double s;
    // location parameter - location of peak - 0 for standard
    double t;
};

/*! 
 * Exponential distribution
 * includes probability density function and cumulative distribution function
 * pdf returns the probability that the variate has the value x
 * cdf returns the range that encompasses x percent of the distribution
 */
class exponential_distribution
{
public:
exponential_distribution(double scale)
        : beta(scale)
    {}
    
    // changing it to assume mu is 0 which is traditionally accepted
    // (1/beta)*(std::exp(-(x-mu)/beta));
    double pdf(double x)
    {
    	return (1/beta)*(std::exp(-x/beta));
    }
    
    // x is the % of the distribution you want covered (e.g for 99% input 0.99)
    // don't take half of this since exponential distribution begins at 0
    double cdf(double x)
    {
    	// std::log is ln
    	if ( beta == 1) {
    		return -std::log(1-x);
    	} else {
    		return -beta * std::log(1-x);
    	}
    }
private:
	// scale parameter - 1 for standard
	double beta;
};

/*! 
 * Dispersal kernel for deterministic spread to cell with highest probability of spread
 *
 * Given a boundaries of the raster, to avoid spreading beyond boundaries of raster
 * Dispersal Kernel type determines use of Exponential or Cauchy distribution
 * to determine probability.
 *
 * Useful for testing as it is deterministic and provides
 * fully replicable results 
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
    
    // TODO could change to store this pair in a tuple
//     std::tuple<int, int> middle;
    int mid_row = 0;
    int mid_col = 0;
    // Can't be in a tuple because it needs to be mutable, but could be in a vector?
    int prev_row = -1;
    int prev_col = -1;
    int row_max;
    int col_max;
    std::vector<std::vector<double>> probability;
    std::vector<std::vector<double>> probability_copy;
    cauchy_distribution cauchy;
    exponential_distribution exponential;
public:
    DeterministicDispersalKernel(DispersalKernelType dispersal_kernel, double distance_scale,
    							double locator, Raster<int> dispersers, double window_size,
    							double ew_res, double ns_res)
        :
        dispersers_(dispersers),
        cauchy(distance_scale, locator),
        exponential(distance_scale),
        kernel_type_(dispersal_kernel),
        east_west_resolution(ew_res),
        north_south_resolution(ns_res)
    {
    	if (kernel_type_ == DispersalKernelType::Cauchy) {
    		// window_size is the percent of all possible dispersal to be included 
    		// in the moving window size
    		row_max = (int) std::ceil(cauchy.cdf(window_size));
    	} else { 
    		// otherwise just default to exponential
    		// if changed this would need to include some check for kernel type
    		row_max = (int) std::ceil(exponential.cdf(window_size));
    	}
    	col_max = row_max / east_west_resolution;
    	row_max /= north_south_resolution;
    	// must be odd
    	if ( row_max % 2 == 0 ) {
    		row_max += 1;
    	}
    	if ( col_max % 2 == 0 ) {
    		col_max += 1;
    	}
    	// create probability matrix
    	// compute distance from center to all cells within window
    	probability.resize(row_max, std::vector<double>(col_max));
    	// also set the size of the copy
    	probability_copy.resize(row_max, std::vector<double>(col_max));
    	mid_row = row_max / 2;
    	mid_col = col_max / 2;
    	double sum = 0;
    	for ( int i = 0; i < row_max; i++ ) {
    		for ( int j = 0; j < col_max; j++ ) {
    			if ( i == mid_row && j == mid_col ) {
    				// do not want anything put in the middle cell
    				probability[i][j] = 0;
    			} else {
    				// store distance to center
    				probability[i][j] = (std::abs(mid_row - i) / east_west_resolution)
    								  + (std::abs(mid_col - j) / north_south_resolution);
    				// determine probability based on distance
    				if (kernel_type_ == DispersalKernelType::Cauchy) {
 					probability[i][j] = std::abs(cauchy.pdf(probability[i][j]));
 					} else if ( kernel_type_ == DispersalKernelType::Exponential ) {
 						probability[i][j] = std::abs(exponential.pdf(probability[i][j]));
 					}
 					sum += probability[i][j];
 				}
    		}
    	}
    	//normalize based on the sum of all probabilities in the raster
    	for ( int i = 0; i < row_max; i++ ) {
    		for ( int j = 0; j < col_max; j++ ) {
    			probability[i][j] /= sum;
    		}
    	}
    }

    /*! Generates a new position for the spread.
     *
     *  Creates smaller matrix "window" from probability matrix
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
    		// determine dimensions for this window
    		
    		// determine the number of dispersers per cell and store these values
    		// in the probability copy
    		// dimensions from the original probability  matrix
    		for ( int i = 0; i < row_max; i++) {
    			for ( int j = 0; j < col_max; j++) {
    				probability_copy[i][j] = probability[i][j];
    			}
    		}
    	}
    	
    	int row_translate = row - mid_row;
    	int col_translate = col - mid_col;
    	
    	double max = (double)-std::numeric_limits<int>::max();
    	// have to add a max value because if row == probability_r_min and the probability
    	// of this cell is greater than every other cell - no cell will be selected
    	int max_prob_row = 0;
    	int max_prob_col = 0;
    	
    	//find cell with highest probability
    	for ( int i = 0; i < row_max; i++ ) {
    		for ( int j = 0; j < col_max; j++ ) {
    			// compare row/col to middle row/col
    			if ((i + row_translate)== row && ((j + col_translate) == col)) {
    				// do nothing
    				// should be zero anyways - could probably remove this
    			}
    			else if (probability_copy[i][j] > max ) {
    				max = probability_copy[i][j];
    				max_prob_row = i;
					max_prob_col = j;
    			}
    		}
    	}

		// subtracting 1 preserves the probability order but ensures that this cell
		// is not select again unless all other cells in window have been selected
    	probability_copy[max_prob_row][max_prob_col] -= 1;
    	
    	prev_row = row;
    	prev_col = col;
    	
    	// need to return values in terms of actual location
    	// find the difference between the row/col values from the middle cell
        return std::make_tuple(max_prob_row + row_translate, max_prob_col + col_translate);
    }

    /*! \copydoc RadialDispersalKernel::supports_kernel()
     */
     // copied from radial_kernel
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
