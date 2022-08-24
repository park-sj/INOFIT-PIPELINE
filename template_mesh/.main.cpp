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

  /* create example Eigen matrix*/
  Eigen::SparseMatrix<double,Eigen::RowMajor> M(4,4);
  std::vector<Eigen::Triplet<double> > tripletList;
  tripletList.push_back(Eigen::Triplet<double>(0,0,1));
  tripletList.push_back(Eigen::Triplet<double>(1,1,1));
  tripletList.push_back(Eigen::Triplet<double>(2,0,1));
  tripletList.push_back(Eigen::Triplet<double>(2,2,1));
  tripletList.push_back(Eigen::Triplet<double>(3,0,1));
  tripletList.push_back(Eigen::Triplet<double>(3,2,1));
  tripletList.push_back(Eigen::Triplet<double>(3,3,1));

  M.setFromTriplets(tripletList.begin(),tripletList.end());

  
  Eigen::SparseMatrix<double>::StorageIndex *outer = M.outerIndexPtr();
  Eigen::SparseMatrix<double>::StorageIndex *inner = M.innerIndexPtr();
  Eigen::SparseMatrix<double>::StorageIndex *innerNonZero = M.innerNonZeroPtr();
  Eigen::SparseMatrix<double>::Scalar *values = M.valuePtr();
  Eigen::Index nonzeros = M.nonZeros();

  std::cout<<"nonzeros"<<std::endl;

  std::cout<<"outer"<<std::endl;
  for (int i=0; i <= M.outerSize(); i++){
    std::cout<<outer[i]<<" ";
  }
  cout<<std::endl;

  std::cout<<"inner"<<std::endl;
  for (int i=0; i <= M.nonZeros() ; i++){
    std::cout<<inner[i]<<" ";
  }
  cout<<std::endl;
 
  std::cout<<"M.innerSize():"<<M.innerSize()<<std::endl;
  std::cout<<"M.nonzeros:"<<M.nonZeros()<<std::endl;
   

  std::cout<<M<<std::endl;
 
  int *outer_arr = (int*)outer;
  int *inner_arr= (int*)inner;
  double *value_arr = (double*)values;

  std::cout<<"outer"<<std::endl;
  for (int i=0; i <= M.outerSize(); i++){
    outer[i]+=1;
    std::cout<<outer[i]<<" ";
  }
  cout<<std::endl;

  std::cout<<"inner"<<std::endl;
  for (int i=0; i <= M.nonZeros() ; i++){
    inner[i]+=1;
    std::cout<<inner[i]<<" ";
  }
  cout<<std::endl;

  std::cout<<"values"<<std::endl;
  for (int i=0; i<7;i++)
  {
    std::cout<<values[i]<<std::endl;
  }
 

  A.num_rows = 4;
  A.num_cols = 4; 
  A.Ap = outer;
  A.Aj = inner;
  
  A.Ax = values;
  A.num_nonzeros = 7;
 
  vec<int, double> b;

  double b_data[5] = {1.0,1.0,1.0,1.0,1.0};
  b.val = b_data;
  b.len = 4;
  
  vec<int,double> x = new_vec<int,double>(4); 

  //gpucg_solve_(*row_ptr, &size_row_ptr, *col_idx, &size_col_idx, *val, &size_val, *rhs, &size_rhs, *x)

  cg<int,double>(A,x,b);

  for(int i =0; i<4; i++)
  {
    std::cout<<x.val[i]<<std::endl;
  }

}
