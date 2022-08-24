#include "sparse_io.h"
#include "vec.h"
#include "gpu_solve.h"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

#define MAX_ERROR 1.0e-6

using namespace std;

template <class IndexType, class ValueType>
void cg(csr_matrix<IndexType, ValueType> A, vec<IndexType, ValueType> x, 
        vec<IndexType, ValueType> b)
{
  int size_b = A.num_rows;
  int size_row_ptr = A.num_rows+1;
  int size_col_idx = A.num_nonzeros;
  gpucg_solve_(A.Ap, &size_row_ptr, A.Aj, &size_col_idx, A.Ax, &size_col_idx, 
               b.val, &size_b, x.val);
}

int main(int argc, char **argv)
{

  csr_matrix<int,double> A;

  /* need to set 
   * num_rows, num_cols, nunm_nonzeros
   * Ap // row pointer
   * Aj // column indices
   * Ax // nonzeros
   */

  int rows[5] = {1,2,3,4,5};
  int cols[4] = {1,2,3,4};

  
  //int rows[5] = {1,2,3,5,8};
  //int cols[7] = {1,2,1,3,1,3,4};

  double values[5] = {1.0,1.0,1.0,1.0,1.0};
 // double values[7] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  A.num_rows = 4;
  A.num_cols = 4; 
  A.num_nonzeros = 4;

  A.Ap = rows;
  A.Aj = cols;
  A.Ax = values;
 
  vec<int, double> b;

  double b_data[4] = {7, 3, 3, 5};
  b.val = b_data;
  b.len = 4;
  
  vec<int,double> x = new_vec<int,double>(4); 
  cg<int,double>(A,x,b);

  for(int i =0; i<4; i++)
  {
    std::cout<<x.val[i]<<std::endl;
  }

}
