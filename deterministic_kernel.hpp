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
 * Cauchy probability density function
 * returns the probability that the variate has the value x
 */
class cauchy_density_function
{
public:
	cauchy_density_function(double scale, double locator)
        : s(scale), t(locator)
    {}

    double operator ()(double x)
    {
    	return 1 / ((s*PI)*(1 + (std::pow((x - t)/s, 2))));
    
    }
private:
	// scale parameter - specifies half width at half maximum - 1 for standard
	// often referred to as lamba which equals 1/beta
    double s;
    // location parameter - location of peak - 0 for standard
    double t;
};

// 
/*! 
 * Exponential probability density function 
 * returns the probability that the variate has the value x
 */
class exponential_density_function
{
public:
exponential_density_function(double scale, double locator)
        : beta(scale), mu(locator)
    {}
    
    double operator ()(double x)
    {
    	return (1/beta)*(std::exp(-(x-mu)/beta));
    }
private:
	// scale parameter - 1 for standard
	double beta;
	// location parameter - 0 for standard 
    double mu;
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
    int row_max_;
    int col_max_;
    Raster<int> dispersers_;
    DispersalKernelType kernel_type_;
    
    int prev_row = -1;
    int prev_col = -1;
    int window_r_min;
    int window_r_max;
    int window_c_min;
    int window_c_max;
    std::vector<std::vector<double>> probability;
    std::vector<std::vector<double>> window;
    cauchy_density_function cauchy_pdf;
    exponential_density_function exponential_pdf;
public:
    DeterministicDispersalKernel(int row_max, int col_max, DispersalKernelType dispersal_kernel,
                          		double distance_scale, double locator, Raster<int> dispersers)
        :
        row_max_(row_max),
        col_max_(col_max),
        dispersers_(dispersers),
        cauchy_pdf(distance_scale, locator),
        exponential_pdf(distance_scale, locator),
        kernel_type_(dispersal_kernel)
    {
    	// create probability matrix
    	// compute distance from center to all cells
    	probability.resize(row_max, std::vector<double>(col_max));
    	int med_row = row_max_/2;
    	int med_col = col_max_/2;
    	double sum = 0;
    	for ( int i = 0; i < row_max_; i++ ) {
    		for ( int j = 0; j < col_max_; j++ ) {
    			// store distance to center
    			probability[i][j] = std::abs(med_row - i) + std::abs(med_col - j);
    			// determine probability based on distance
    			if (kernel_type_ == DispersalKernelType::Cauchy) {
 					probability[i][j] = std::abs(cauchy_pdf(probability[i][j]));
 				} else if ( kernel_type_ == DispersalKernelType::Exponential ) {
 					probability[i][j] = std::abs(exponential_pdf(probability[i][j]));
 				}
 				sum += probability[i][j];
    		}
    	}
    	//normalize based on the sum of all probabilities in the raster
    	for ( int i = 0; i < row_max_; i++ ) {
    		for ( int j = 0; j < col_max_; j++ ) {
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
    		int dispersers_num = dispersers_(row, col)/2;
    		
    		// determine moving window boundary values
    		// length/width = number of dispersers
    		// TODO could change to different dimensions
    		window_r_min = ((row - (dispersers_num)) < 0) ? 0 : (row - (dispersers_num));
    		window_r_max = ((row + (dispersers_num)) > row_max_) ? row_max_ : (row + (dispersers_num));
			window_c_min = ((col - (dispersers_num)) < 0 ) ? 0 : (col - (dispersers_num));
    		window_c_max = ((col + (dispersers_num)) > col_max_) ? col_max_ : (col + (dispersers_num));
    		window.resize((window_r_max - window_r_min), std::vector<double>(window_c_max - window_c_min));
    		// copy probabilities from window_r_max matrix
    		for ( int i = window_r_min; i < window_r_max; i++) {
    			for ( int j = window_c_min; j < window_c_max; j++) {
    				window[i][j] = probability[i][j];
    			}
    		}
    	}
    	
    	double max = (double)-std::numeric_limits<int>::max();
    	// have to add a max value because if row == window_r_min and the probability
    	// of this cell is greater than every other cell - no cell will be selected
    	int max_prob_row = window_r_min;
    	int max_prob_col = window_c_min;
    	
    	//find cell with highest probability
    	for ( int i = window_r_min; i < window_r_max; i++ ) {
    		for ( int j = window_c_min; j < window_c_max; j++ ) {
    			if (i == row && j == col) {
    				// do nothing
    			}
    			else if (window[i][j] > max ) {
    				max = window[i][j];
    				max_prob_row = i;
					max_prob_col = j;
    			}
    		}
    	}

		// subtracting 1 preserves the probability order but ensures that this cell
		// is not select again unless all other cells in window have been selected
    	window[max_prob_row][max_prob_col] -= 1;
    	
    	prev_row = row;
    	prev_col = col;
        return std::make_tuple(max_prob_row, max_prob_col);
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
