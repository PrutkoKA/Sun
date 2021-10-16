#pragma once

#include <vector>
#include <iostream>

class mass_matrix
{
public:
	std::vector < double > a_diagonal;
	std::vector < double > b_diagonal;
	std::vector < double > c_diagonal;
	std::vector < std::vector < double > > inversed_mass_matrix;

	void calculate_mass_matrix(const std::vector < double > &x);
	void calculate_mass_matrix(std::vector < double >::const_iterator vec_begin, std::vector < double >::const_iterator vec_end, unsigned int vec_size);
	void fill_inverse_mass_matrix();
	void print_inverse_matrix();
};