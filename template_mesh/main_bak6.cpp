#include <iostream>
#include <fstream>
#include <vector>
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
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Iterator_range.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SVD>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

#include "sparse_io.h"
#include "vec.h"
#include "gpu_solve.h"

#include <time.h>

#include <string.h>
#include <stdio.h>



typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Segment_3 Segment;
typedef K::Ray_3 Ray;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;

/* functions to define 
1. Create D, 
 2. Create U.
 3. Create M
*/


template<typename T,typename S>
void create_incidence_matrix(CGAL::Surface_mesh<T>& Mesh, Eigen::SparseMatrix<S, Eigen::RowMajor> &M, double alpha); 

template<typename T>
void stack_sparse_vertical(Eigen::SparseMatrix<T, Eigen::RowMajor>& upper, Eigen::SparseMatrix<T, Eigen::RowMajor>& lower, Eigen::SparseMatrix<T, Eigen::RowMajor>& stacked );

template<typename T,typename S,typename R>
int create_correspondence_matrix(CGAL::AABB_tree<R>& tree,CGAL::Surface_mesh<T>& Mesh, Eigen::SparseMatrix<S, Eigen::RowMajor>& D, Eigen::SparseMatrix<S, Eigen::RowMajor>& U);

template<typename T, typename S>
void allocate_cuda_matrix(const  Eigen::SparseMatrix<S, Eigen::RowMajor>& A_in, csr_matrix<T, S>& A_out);

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

int main(int argc, char *argv[])
{   
	//Pre
    Mesh templateMesh, targetMesh;

    std::cout<<argv[1] << "\n";
    std::ifstream templateInput(argv[1]);

    std::cout<<argv[2] << "\n";
    std::ifstream targetInput(argv[2]);

    std::cout<<argv[3] << "\n";
    std::ofstream pointsout(argv[3]);
    
    templateInput >> templateMesh;
    targetInput >> targetMesh;

    std::vector< std::tuple<size_t, size_t> > directedEdges;



      /* construct tree on target matrix */
      Mesh::Face_range face_r = targetMesh.faces();
      Mesh::Face_range::iterator f_begin, f_end, f_itr;
      f_begin = face_r.begin();
      f_end = face_r.end(); 

      Tree targetMeshTree(faces(targetMesh).first, faces(targetMesh).second, targetMesh);
      try
      {
	  targetMeshTree.accelerate_distance_queries();  

      } 
      catch (int e)
      {
	  std::cout<<e<<std::endl;
	  return -1;
      }    

      /* concludes preprocessing */
      
      Eigen::SparseMatrix<double, Eigen::RowMajor> M;
      create_incidence_matrix<Point, double>(templateMesh, M, 1000);
      
      /* preallocate D,U */
	size_t template_vertex_count = templateMesh.number_of_vertices();

  Eigen::SparseMatrix<double, Eigen::RowMajor> somematrix;

	Eigen::SparseMatrix<double, Eigen::RowMajor> D(template_vertex_count, 4* template_vertex_count);
	Eigen::SparseMatrix<double, Eigen::RowMajor> U(template_vertex_count, 3);
  
      /* also allocater */
	create_correspondence_matrix<Point, double,Traits>(targetMeshTree, templateMesh, D,U);
	
	/* Concatenate M and D */
	Eigen::SparseMatrix<double, Eigen::RowMajor> A;
	stack_sparse_vertical<double>(M,D,A);

	/* Allocate corresponding Zero matrix for M */
	size_t M_rows = M.rows();
	Eigen::SparseMatrix<double, Eigen::RowMajor> Z(M_rows,3);
	
	Eigen::SparseMatrix<double, Eigen::RowMajor> B;
	stack_sparse_vertical<double>(Z,U, B);

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
	
	Eigen::SparseMatrix<double, Eigen::RowMajor> AtA = Eigen::SparseMatrix<double, Eigen::RowMajor>(A.transpose()) * A;
	//Eigen::MatrixX3D AtB = Eigen::SparseMatrix<double, Eigen::RowMajor>(A.transpose())*B;
	Eigen::SparseMatrix<double, Eigen::RowMajor> AtB = Eigen::SparseMatrix<double, Eigen::RowMajor>(A.transpose())*B;
	
	std::cout<<"solving"<<std::endl;

	std::cout<<"AtA.rows():"<<AtA.rows()<<std::endl;
		std::cout<<"AtA.cols():"<<AtA.cols()<<std::endl;
	
	AtA.makeCompressed(); 
	//solver.compute(AtA);
	//Eigen::MatrixXd X = solver.solve(AtB);

	Eigen::MatrixXd col_1 = AtB.col(0);
  Eigen::MatrixXd col_2 = AtB.col(1);
  Eigen::MatrixXd col_3 = AtB.col(2); 

	double *col_1_values = col_1.data();
  double *col_2_values = col_2.data();
  double *col_3_values = col_3.data(); 

	int col_length = AtB.rows();

	vec<int, double> b_1;
  vec<int, double> b_2;
  vec<int, double> b_3;

	b_1.len = col_length;
	b_1.val = col_1_values;

  b_2.len = col_length;
  b_2.val = col_2_values;

  b_3.len = col_length;
  b_3.val = col_3_values;

    int outerStartsLength = AtA.outerSize();
    int valuesLength = AtA.nonZeros();

  std::cout<<outerStartsLength<<std::endl;
  std::cout<<valuesLength<<std::endl;

	csr_matrix<int, double> AtA_cuda_1;
  allocate_cuda_matrix<int, double>(AtA, AtA_cuda_1);
	vec<int, double> x_1 = new_vec<int,double>(col_length);
	cg<int, double>(AtA_cuda_1, x_1, b_1);

	std::cout<<x_1.val[1]<<std::endl;
	std::cout<<x_1.val[2]<<std::endl;

  std::cout<<"**************"<<std::endl;
  AtA = Eigen::SparseMatrix<double, Eigen::RowMajor>(A.transpose()) * A;
  outerStartsLength = AtA.outerSize();
  valuesLength = AtA.nonZeros();

  std::cout<<outerStartsLength<<std::endl;
  std::cout<<valuesLength<<std::endl;

  csr_matrix<int, double> AtA_cuda_2;
  allocate_cuda_matrix<int, double>(AtA, AtA_cuda_2);
  vec<int,double> x_2 = new_vec<int, double>(col_length);
  cg<int,double>(AtA_cuda_2, x_2, b_2);

  std::cout<<x_2.val[1]<<std::endl;
  std::cout<<x_2.val[2]<<std::endl;

  std::cout<<"**************"<<std::endl;
  AtA = Eigen::SparseMatrix<double, Eigen::RowMajor>(A.transpose()) * A;
  outerStartsLength = AtA.outerSize();
  valuesLength = AtA.nonZeros();

  std::cout<<outerStartsLength<<std::endl;
  std::cout<<valuesLength<<std::endl;

  csr_matrix<int, double> AtA_cuda_3;
  allocate_cuda_matrix<int, double>(AtA, AtA_cuda_3);
  vec<int,double> x_3 = new_vec<int, double>(col_length);
  cg<int,double>(AtA_cuda_3, x_3, b_3);

  std::cout<<x_3.val[1]<<std::endl;
  std::cout<<x_3.val[2]<<std::endl;

  /* create matrix X from x_1,x_2,x_3 */

  

	/* iterate through mesh and update    
	Mesh::Vertex_range vertex_r = templateMesh.vertices();
	Mesh::Vertex_range::iterator v_begin, v_end, v_itr; 

	v_begin = vertex_r.begin();
	v_end = vertex_r.end();
	
	
	for( v_itr = v_begin; v_itr != v_end; v_itr++ )
	{
	    Point p = templateMesh.point(*v_itr);
	    size_t idx = *v_itr;
	    std::cout<<idx<<std::endl;
	    double v1 = p[0]; 
	    double v2 = p[1];
	    double v3 = p[2]; 
	    double v4 = 1.0;
	    Eigen::Matrix<double,1,4> v;
	    v << v1,v2,v3,v4;
	    //Eigen::Vector4f v(v1,v2,v3,v4);
	    Eigen::MatrixXd Xsub = X.block<4,3>(4*idx,0);
	    Eigen::ArrayXd Af;
	    Af = Eigen::VectorXd(U.row(idx));
	    Eigen::VectorXd vt= v.transpose();
	    Eigen::VectorXd u_approx = v*Xsub;
	    Point p2(u_approx[0], u_approx[1], u_approx[2]);
	    templateMesh.point(*v_itr) = p2; 
	}*/
	
    pointsout<<templateMesh;

    return EXIT_SUCCESS;
}

template<typename T, typename S>
void allocate_cuda_matrix(const Eigen::SparseMatrix<S, Eigen::RowMajor>& A_in, csr_matrix<T, S>& A_out)
{
    const Eigen::SparseMatrix<double, Eigen::RowMajor>::StorageIndex *outerStarts = A_in.outerIndexPtr();
    const Eigen::SparseMatrix<double, Eigen::RowMajor>::StorageIndex *innerIndices = A_in.innerIndexPtr();
    const Eigen::SparseMatrix<double, Eigen::RowMajor>::Scalar *values = A_in.valuePtr();


    /* copy three arrays into new memory */
    int outerStartsLength = A_in.outerSize();
    int valuesLength = A_in.nonZeros();
  
    int* csr_row_pointers = new int[outerStartsLength];
    int* csr_col_pointers = new int[valuesLength]; 
    double* csr_values = new double[valuesLength];

    /* add one to index to convert to one based indexing */
    for (int i = 0; i <A_in.outerSize(); i++)
    {
        csr_row_pointers[i] = outerStarts[i]+1;

    }
    for (int i = 0; i <A_in.nonZeros(); i++)
    {
        csr_col_pointers[i] = innerIndices[i]+1;
        csr_values[i] = values[i];
    }


    A_out.num_rows =(int)A_in.rows();
    A_out.num_cols =(int)A_in.cols();
    A_out.Ax = csr_values;
    A_out.num_nonzeros = valuesLength;
    A_out.Ap = csr_row_pointers;
    A_out.Aj = csr_col_pointers;
}


template<typename T,typename S>
void create_incidence_matrix(CGAL::Surface_mesh<T>& mesh, Eigen::SparseMatrix<S, Eigen::RowMajor>& out, double alpha)
{
    typedef typename::CGAL::Surface_mesh<T> Mesh;
    typename::CGAL::Surface_mesh<T>::Face_range mesh_face_r = mesh.faces();
    typename::CGAL::Surface_mesh<T>::Face_range::iterator f_begin, f_end, f_itr;
    std::vector< std::tuple<size_t, size_t> > directedEdges;
    std::vector< Eigen::Triplet<double> > MtripletList;

    f_begin = mesh_face_r.begin();
    f_end = mesh_face_r.end();
     
    for( f_itr = f_begin; f_itr != f_end; f_itr++ )
    { 
        typename::CGAL::Surface_mesh<T>::Face_index f_dr = *f_itr;
        CGAL::Iterator_range<CGAL::Halfedge_around_face_iterator<Mesh>> face_halfedge_r
		= halfedges_around_face( mesh.halfedge(f_dr), mesh);

        CGAL::Halfedge_around_face_iterator<Mesh> e_begin, e_end, e_itr;
        e_begin = face_halfedge_r.begin();
        e_end = face_halfedge_r.end();       

        for( e_itr = e_begin; e_itr != e_end; e_itr++ )
        {
            halfedge_descriptor e_descriptor = *e_itr; 
            typename::Mesh::Edge_range edges = mesh.edges();
      
            size_t src_idx = source(e_descriptor, mesh);
            size_t dst_idx = target(e_descriptor, mesh);

            if (src_idx < dst_idx ){
                directedEdges.push_back(std::tuple<size_t, size_t>(src_idx, dst_idx));    
            }

        }
    }
    std::cout<<directedEdges.size()<<std::endl;

    std::vector< std::tuple<size_t, size_t> >::iterator itr;
    itr =directedEdges.begin();
    size_t row_num = 0;

     for (itr ; itr != directedEdges.end(); itr++)
     {
        std::tuple<size_t, size_t> edge = *itr;
        size_t src_idx = get<0>(edge);
        size_t dst_idx = get<1>(edge);

        for(int i = 0; i<4; i++)
        {
            int j = 4*row_num+i;
            int k = 4*src_idx+i; 
            int l = 4*dst_idx+i; 
            MtripletList.push_back(Eigen::Triplet<double>(j,k,-1*alpha));
            MtripletList.push_back(Eigen::Triplet<double>(j,l,1*alpha));
        }
        row_num++; 
    }
    size_t N, M; 
    N = directedEdges.size(); // row number 
    M = mesh.number_of_vertices();    // vertex number 
    Eigen::SparseMatrix<S, Eigen::RowMajor> incidenceMatrix(4*N,4*M);
    incidenceMatrix.setFromTriplets(MtripletList.begin(), MtripletList.end()); 
    out = incidenceMatrix;
}

template<typename T,typename S, typename R>
int create_correspondence_matrix(CGAL::AABB_tree<R>& tree, CGAL::Surface_mesh<T>& mesh, Eigen::SparseMatrix<S, Eigen::RowMajor>& D, Eigen::SparseMatrix<S, Eigen::RowMajor>& U)
{
   /* Iterate through template mesh and find closest poin on target mesh*/
    typedef typename::CGAL::Surface_mesh<T> Mesh;
    typename::Mesh::Vertex_range vertex_r = mesh.vertices();
    typename::Mesh::Vertex_range::iterator v_begin, v_end, v_itr; 
    v_begin = vertex_r.begin();
    v_end = vertex_r.end();
    std::vector< Eigen::Triplet<double> > DtripletList;
    std::vector< Eigen::Triplet<double> > UtripletList; 

    /*
    std::vector<double> Ucol1;
    std::vector<double> Ucol2;
    std::vector<double> Ucol3;
    */
    size_t row = 0;

    for( v_itr = v_begin; v_itr != v_end; v_itr++ )
    {
	int v1_idx, v2_idx, v3_idx, v4_idx;
	
	size_t v_idx = *v_itr; 
	
	/* if v_idx is one of the landmark indices, process separately*/

	
	Point p = mesh.point(*v_itr);
		
	Point_and_primitive_id pp = tree.closest_point_and_primitive(p);
	Point closest = pp.first;
	v1_idx = 4*row;
	v2_idx = 4*row+1;
	v3_idx = 4*row+2;
	v4_idx = 4*row+3;
	 
	double v1 = p[0]; 
  double v2 = p[1];
  double v3 = p[2]; 
  double v4 = 1.0;

  double u1 = closest[0];
  double u2 = closest[1];
  double u3 = closest[2];
           
  DtripletList.push_back(Eigen::Triplet<double>(row, v1_idx, v1));
  DtripletList.push_back(Eigen::Triplet<double>(row, v2_idx, v2));
  DtripletList.push_back(Eigen::Triplet<double>(row, v3_idx, v3));
  DtripletList.push_back(Eigen::Triplet<double>(row, v4_idx, v4));

  UtripletList.push_back(Eigen::Triplet<double>(row, 0, u1));
  UtripletList.push_back(Eigen::Triplet<double>(row, 1, u2));
  UtripletList.push_back(Eigen::Triplet<double>(row, 2, u3));
  
	/*
        Ucol1.push_back(u1);
        Ucol2.push_back(u2);
        Ucol3.push_back(u3);
	*/
        row++;
    }
  D.setFromTriplets(DtripletList.begin(), DtripletList.end());
  U.setFromTriplets(UtripletList.begin(), UtripletList.end());
    
    return 0;
}



template<typename T>
void stack_sparse_vertical(Eigen::SparseMatrix<T, Eigen::RowMajor>& upper, Eigen::SparseMatrix<T, Eigen::RowMajor>& lower, Eigen::SparseMatrix<T, Eigen::RowMajor>& stacked)
{
    std::vector< Eigen::Triplet<T> > tripletList;

    for( int k = 0; k < upper.outerSize(); ++k)
    {
      for(typename::Eigen::SparseMatrix<T, Eigen::RowMajor>::InnerIterator it(upper, k); it; ++it)
      {
        it.value();
        size_t row_idx = it.row(); 
        size_t col_idx = it.col();
        size_t it_idx = it.index();
        double it_value = it.value();
        // fill triplet with this 
        tripletList.push_back(Eigen::Triplet<T>(row_idx, col_idx, it_value));
      }
    }
    
    for( int k = 0; k < lower.outerSize(); ++k)
    {
       for(typename::Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(lower, k); it; ++it)
       {
          it.value();
          size_t row_idx = it.row();
          size_t col_idx = it.col();
          size_t it_idx = it.index();
          double it_value = it.value();
          tripletList.push_back(Eigen::Triplet<double>(row_idx, col_idx, it_value));
       }
     }
  
    size_t stacked_rows = upper.rows()+lower.rows();
    size_t stacked_cols = upper.cols();
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(stacked_rows, stacked_cols);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    stacked=A;
}
