#include "mass_matrix.h"

void mass_matrix::calculate_mass_matrix(const std::vector < double >& x)
{
	calculate_mass_matrix(x.begin(), x.end(), x.size());
}

void mass_matrix::calculate_mass_matrix(std::vector < double >::const_iterator vec_begin, std::vector < double >::const_iterator vec_end, unsigned int vec_size)
{
	a_diagonal.resize(vec_size);
	b_diagonal.resize(vec_size);
	c_diagonal.resize(vec_size);
	a_diagonal[vec_size - 1] = 0.25;
	c_diagonal[0] = 0.25;
	unsigned int i = 0;
	std::vector < double >::const_iterator mid_it, right_it;
	for (auto it = vec_begin; it != vec_end; ++it)
	{
		auto left_it = mid_it = right_it = it;
		if (it != vec_begin)
			--left_it;
		
		++right_it;

		b_diagonal[i] = 0.75;
		if (it != vec_begin && right_it != vec_end)
		{
			a_diagonal[i] = 0.25 * (*mid_it - *left_it) / (*right_it - *left_it);
			c_diagonal[i] = 0.25 * (*right_it - *mid_it) / (*right_it - *left_it);
		}
		++i;
	}
}

void mass_matrix::fill_inverse_mass_matrix()
{
	std::vector < double > a(a_diagonal), b(b_diagonal), c(c_diagonal);
	unsigned int n = a_diagonal.size() - 1;

	inversed_mass_matrix.resize(a_diagonal.size());
	for (int i = 0; i < a_diagonal.size(); ++i)
		inversed_mass_matrix[i].resize(a_diagonal.size());

	for (int inv_m_col = 0; inv_m_col < a_diagonal.size(); ++inv_m_col)
	{
		std::vector < double > a(a_diagonal), b(b_diagonal), c(c_diagonal);
		std::vector< double > d(a_diagonal.size());

		d[inv_m_col] = 1.;

		c[0] /= b[0];
		d[0] /= b[0];

		for (int i = 1; i < n; ++i)
		{
			c[i] /= b[i] - a[i] * c[i - 1];
			d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
		}

		d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);
		inversed_mass_matrix[n][inv_m_col] = d[n];

		for (int i = n; i-- > 0;) 
		{
			d[i] -= c[i] * d[i + 1];
			inversed_mass_matrix[i][inv_m_col] = d[i];		// i is row for inversed mass matrix
		}
	}
}

void mass_matrix::print_inverse_matrix()
{
	for (int i = 0; i < inversed_mass_matrix.size(); ++i)
	{
		for (int j = 0; j < inversed_mass_matrix[i].size(); ++j)
		{
			std::cout << inversed_mass_matrix[i][j] << "\t";
		}
		std::cout << "\n";
	}
}