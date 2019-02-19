#ifndef FUSION_NODE__HUNGARIAN_HPP_
#define FUSION_NODE__HUNGARIAN_HPP_

#include <iostream>
#include <vector>

using namespace std;


class HungarianAlgorithm
{
public:
	HungarianAlgorithm();
	~HungarianAlgorithm();
	void solve(vector<vector<double> >& dist_matrix, vector<int>& out_assignment,
               unsigned int n_rows, unsigned int n_cols);

private:
  int m_n_rows;
  int m_n_cols;
  int m_n_elements;
  double m_dist_matrix[256 * 256] = {0};
  double m_assignment[256] = {0};
  bool m_covered_columns[256] = {};
  bool m_covered_rows[256] = {};
  bool m_star_matrix[256 * 256] = {};
  bool m_prime_matrix[256 * 256] = {};
  bool m_new_star_matrix[256 * 256] = {};

  void set_col_row(int n_rows, int n_cols) {
    m_n_rows = n_rows;
    m_n_cols = n_cols;
    m_n_elements = n_rows * n_cols;
  }
  void assignmentoptimal();
	void buildassignmentvector();
	void step2a(int minDim);
	void step2b(int minDim);
	void step3(int minDim);
	void step4(int minDim, int row, int col);
	void step5(int minDim);
};

#endif // FUSION_NODE__FUSION_NODE_HPP_

