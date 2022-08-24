#include "sparse_io.h"
#include "vec.h"
#include "gpu_solve.h"
#include <iostream>


#include <iostream>
#include <fstream>
#include <tuple>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Surface_mesh.h>
#include <string>
#include <vector>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
//#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Iterator_range.h>
//#include <CGAL/Eigen_matrix.h>
//#include <CGAL/KroneckerProduct>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SVD>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

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

 Eigen::SparseMatrix<double, Eigen::RowMajor> M(2,2);
 std::vector<Eigen::Triplet<double> > tripletList;
 tripletList.push_back(Eigen::Triplet<double>(0,0,4));
 tripletList.push_back(Eigen::Triplet<double>(0,1,1));
 tripletList.push_back(Eigen::Triplet<double>(1,0,1));
 tripletList.push_back(Eigen::Triplet<double>(1,1,3));
 
 M.setFromTriplets(tripletList.begin(), tripletList.end());
 
 Eigen::SparseMatrix<double>::StorageIndex *rows = M.outerIndexPtr();
 Eigen::SparseMatrix<double>::StorageIndex *cols = M.innerIndexPtr();
 Eigen::SparseMatrix<double>::Scalar *values = M.valuePtr();

 std::cout<<"rows:"<<std::endl; 
 for(int i = 0; i <=M.outerSize(); i++)
 {
   rows[i]+=1;
   std::cout<<rows[i]<<" ";
 }
 std::cout<<endl;

 std::cout<<"cols:"<<std::endl;
 for(int i = 0; i <=M.nonZeros(); i++)
 {
   cols[i]+=1;
   std::cout<<cols[i]<<" ";
 }
 std::cout<<endl;

 std::cout<<"values"<<std::endl;
 for(int i = 0; i <=M.nonZeros(); i++)
 {
    std::cout<<values[i]<<" ";
 }
 std::cout<<endl;

 
 std::cout<<M<<std::endl;

 /*
 double values[4] = {4,1,1,3};
 int rows[3] = {1,3,5};
 int cols[4] = {1,2,1,2};
*/
  A.num_rows = 2;
  A.num_cols = 2; 
  A.num_nonzeros = 4;

  A.Ap = rows;
  A.Aj = cols;
  
  A.Ax = values;
 
  vec<int, double> b;

  double b_data[2] = {1, 2};
  b.val = b_data;
  b.len = 2;
  
  vec<int,double> x = new_vec<int,double>(2); 
  cg<int,double>(A,x,b);

  for(int i =0; i<2; i++)
  {
    std::cout<<x.val[i]<<std::endl;
  }

}
