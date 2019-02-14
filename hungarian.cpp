#include <stdlib.h>
#include <cfloat> // for DBL_MAX
#include <cmath>  // for fabs()
#include "hungarian.hpp"


HungarianAlgorithm::HungarianAlgorithm(): m_n_rows(),
                                          m_n_cols()
{}

HungarianAlgorithm::~HungarianAlgorithm(){}


//********************************************************//
// A single function wrapper for solving assignment problem.
//********************************************************//
void HungarianAlgorithm::solve(vector<vector<double> >& dist_matrix,
                               vector<int>& out_assignment,
                               unsigned int n_rows, unsigned int n_cols)
{
	// Fill in the m_dist_matrix. Mind the index is "i + m_n_rows * j".
	// Here the cost matrix of size MxN is defined as a double precision array of N*M elements.
	// In the solving functions matrices are seen to be saved MATLAB-internally in row-order.
	// (i.e. the matrix [1 2; 3 4] will be stored as a vector [1 3 2 4], NOT [1 2 3 4]).
	for (int i = 0; i < n_rows; ++i) {
		for (int j = 0; j < n_cols; ++j) {
      double value = dist_matrix[i][j];
  		if (value < 0) {
  			cerr << "All matrix elements have to be non-negative." << endl;
      }
			m_dist_matrix[i + n_rows * j] = value;
    }
  }

	// call solving function
  set_col_row(n_rows, n_cols);
	assignmentoptimal();

	for (int r = 0; r < m_n_rows; ++r) {
    out_assignment[r] = m_assignment[r];
    m_assignment[r] = -1;
  }

  // clean
  for (int i = 0; i < m_n_cols; ++i) {
    m_covered_columns[i] = false;
  }
  for (int i = 0; i < m_n_rows; ++i) {
    m_covered_rows[i] = false;
  }
  for (int i = 0; i < m_n_elements; ++i) {
    m_star_matrix[i] = false;
    m_prime_matrix[i] = false;
    m_new_star_matrix[i] = false;
    m_dist_matrix[i] = 0;
  }
}


//********************************************************//
// Solve optimal solution for assignment problem using Munkres algorithm,
// also known as Hungarian Algorithm.
//********************************************************//
void HungarianAlgorithm::assignmentoptimal()
{
	/* initialization */
  for (int row = 0; row < m_n_rows; ++row) {
		m_assignment[row] = -1;
  }

	/* generate working copy of distance Matrix */
	/* check if all matrix elements are positive */
	double *dist_matrix_end = m_dist_matrix + m_n_elements;

	/* preliminary steps */
  int minDim;
	if (m_n_rows <= m_n_cols) {
		minDim = m_n_rows;

		for (int row = 0; row < m_n_rows; ++row) {
			/* find the smallest element in the row */
			double *tmp_dist_matrix = m_dist_matrix + row;
		  double minValue = *tmp_dist_matrix;
			tmp_dist_matrix += m_n_rows;
			while (tmp_dist_matrix < dist_matrix_end) {
				double value = *tmp_dist_matrix;
				if (value < minValue) {
					minValue = value;
        }
				tmp_dist_matrix += m_n_rows;
			}

			/* subtract the smallest element from each element of the row */
			tmp_dist_matrix = m_dist_matrix + row;
			while (tmp_dist_matrix < dist_matrix_end) {
				*tmp_dist_matrix -= minValue;
				tmp_dist_matrix += m_n_rows;
			}
		}

		/* Steps 1 and 2a */
		for (int row = 0; row < m_n_rows; ++row) {
			for (int col = 0; col < m_n_cols; ++col) {
				if (fabs(m_dist_matrix[row + m_n_rows * col]) < DBL_EPSILON) {
					if (!m_covered_columns[col]) {
						m_star_matrix[row + m_n_rows * col] = true;
						m_covered_columns[col] = true;
						break;
					}
        }
      }
    }
	}	else {
		minDim = m_n_cols;

		for (int col = 0; col < m_n_cols; ++col) {
			/* find the smallest element in the column */
			double *tmp_dist_matrix = m_dist_matrix + m_n_rows * col;
			double *column_end = tmp_dist_matrix + m_n_rows;

			double minValue = *tmp_dist_matrix++;
			while (tmp_dist_matrix < column_end) {
				double value = *tmp_dist_matrix++;
				if (value < minValue) {
					minValue = value;
        }
			}

			/* subtract the smallest element from each element of the column */
			tmp_dist_matrix = m_dist_matrix + m_n_rows * col;
			while (tmp_dist_matrix < column_end) {
				*tmp_dist_matrix++ -= minValue;
      }
		}

		/* Steps 1 and 2a */
		for (int col = 0; col < m_n_cols; ++col) {
			for (int row = 0; row < m_n_rows; ++row) {
				if (fabs(m_dist_matrix[row + m_n_rows * col]) < DBL_EPSILON) {
					if (!m_covered_rows[row]) {
						m_star_matrix[row + m_n_rows * col] = true;
						m_covered_columns[col] = true;
						m_covered_rows[row] = true;
						break;
					}
        }
      }
    }
		for (int row = 0; row < m_n_rows; ++row) {
			m_covered_rows[row] = false;
    }
	}

	/* move to step 2b */
	step2b(minDim);

	return;
}

/********************************************************/
void HungarianAlgorithm::buildassignmentvector()
{
	for (int row = 0; row < m_n_rows; ++row) {
		for (int col = 0; col < m_n_cols; ++col) {
			if (m_star_matrix[row + m_n_rows * col]) {
				m_assignment[row] = col;
				break;
			}
    }
  }
}

/********************************************************/
void HungarianAlgorithm::step2a(int minDim)
{
	/* cover every column containing a starred zero */
	for (int col = 0; col < m_n_cols; ++col) {
		bool *tmp_start_matrix = m_star_matrix + m_n_rows * col;
		bool *column_end = tmp_start_matrix + m_n_rows;
		while (tmp_start_matrix < column_end){
			if (*tmp_start_matrix++) {
				m_covered_columns[col] = true;
				break;
			}
		}
	}
	/* move to step 3 */
	step2b(minDim);
}

/********************************************************/
void HungarianAlgorithm::step2b(int minDim)
{
	/* count covered columns */
	int n_covered_columns = 0;
	for (int col = 0; col < m_n_cols; ++col) {
		if (m_covered_columns[col]) {
			++n_covered_columns;
    }
  }

	if (n_covered_columns == minDim) {
		/* algorithm finished */
		buildassignmentvector();
	} else {
		/* move to step 3 */
		step3(minDim);
	}
}

/********************************************************/
void HungarianAlgorithm::step3(int minDim)
{
	bool zerosFound = true;
	while (zerosFound)
	{
		zerosFound = false;
		for (int col = 0; col < m_n_cols; ++col) {
			if (!m_covered_columns[col]) {
				for (int row = 0; row < m_n_rows; ++row) {
					if ((!m_covered_rows[row]) &&
              (fabs(m_dist_matrix[row + m_n_rows * col]) < DBL_EPSILON)) {
						/* prime zero */
						m_prime_matrix[row + m_n_rows * col] = true;

						/* find starred zero in current row */
            int star_col = 0;
						for (; star_col < m_n_cols; ++star_col) {
							if (m_star_matrix[row + m_n_rows * star_col]) {
								break;
              }
            }

						if (star_col == m_n_cols) { /* no starred zero found */
							/* move to step 4 */
							step4(minDim, row, col);
							return;
						} else {
							m_covered_rows[row] = true;
							m_covered_columns[star_col] = false;
							zerosFound = true;
							break;
						}
					}
        }
      }
    }
	}

	/* move to step 5 */
	step5(minDim);
}

/********************************************************/
void HungarianAlgorithm::step4(int minDim, int row, int col)
{
	/* generate temporary copy of starMatrix */
	for (int n = 0; n < m_n_elements; ++n) {
		m_new_star_matrix[n] = m_star_matrix[n];
  }

	/* star current zero */
	m_new_star_matrix[row + m_n_rows * col] = true;

	/* find starred zero in current column */
	int star_col = col;
  int star_row = 0;
	for (; star_row<m_n_rows; ++star_row) {
		if (m_star_matrix[star_row + m_n_rows * star_col]) {
			break;
    }
  }

	while (star_row < m_n_rows) {
		/* unstar the starred zero */
		m_new_star_matrix[star_row + m_n_rows * star_col] = false;

		/* find primed zero in current row */
		int prime_row = star_row;
    int prime_col = 0;
		for (; prime_col < m_n_cols; ++prime_col) {
			if (m_prime_matrix[prime_row + m_n_rows * prime_col]) {
				break;
      }
    }

		/* star the primed zero */
		m_new_star_matrix[prime_row + m_n_rows * prime_col] = true;

		/* find starred zero in current column */
		star_col = prime_col;
		for (star_row = 0; star_row < m_n_rows; ++star_row) {
			if (m_star_matrix[star_row + m_n_rows * star_col]) {
				break;
      }
    }
	}

	/* use temporary copy as new m_star_matrix */
	/* delete all primes, uncover all rows */
	for (int n = 0; n < m_n_elements; ++n) {
		m_prime_matrix[n] = false;
		m_star_matrix[n] = m_new_star_matrix[n];
	}

	for (int n = 0; n < m_n_rows; ++n) {
		m_covered_rows[n] = false;
  }

	/* move to step 2a */
	step2a(minDim);
}

/********************************************************/
void HungarianAlgorithm::step5(int minDim)
{
	/* find smallest uncovered element h */
	double h = DBL_MAX;
	for (int row = 0; row < m_n_rows; ++row) {
		if (!m_covered_rows[row]) {
			for (int col = 0; col < m_n_cols; ++col) {
				if (!m_covered_columns[col]) {
					double value = m_dist_matrix[row + m_n_rows * col];
					if (value < h) {
						h = value;
          }
				}
      }
    }
  }

	/* add h to each covered row */
	for (int row = 0; row < m_n_rows; ++row) {
		if (m_covered_rows[row]) {
			for (int col = 0; col < m_n_cols; ++col) {
				m_dist_matrix[row + m_n_rows * col] += h;
      }
    }
  }

	/* subtract h from each uncovered column */
	for (int col = 0; col < m_n_cols; ++col) {
		if (!m_covered_columns[col]) {
			for (int row = 0; row < m_n_rows; ++row) {
				m_dist_matrix[row + m_n_rows * col] -= h;
      }
    }
  }

	/* move to step 3 */
	step3(minDim);
}

