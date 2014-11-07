#include "CSR.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_cblas.h>



class Dense_matrix 
{
public:
	gsl_matrix *data_;
	size_t nrow_;
	size_t ncol_;


	// Member functions

	/// Constructor
	Dense_matrix ();
	Dense_matrix (size_t nrow, size_t ncol);
	void init_alloc(size_t nrow, size_t ncol);
	/// Destructor
	~Dense_matrix ();
	/// Initialize the matrix with a particular value
	void initialize_matrix (double val);	
	// 2nd order forbious norm
	float get_frobenious_norm ();

	// int gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C)
	void matrix_multiply (gsl_matrix* b, gsl_matrix* result);

	/// gsl_matrix invert
	/// gsl_permutation * gsl_permutation_alloc (size_t n)
	/// void gsl_permutation_free (gsl_permutation * p)
	/// int gsl_linalg_LU_decomp (gsl_matrix * A, gsl_permutation * p, int * signum)
	/// int gsl_linalg_LU_invert (const gsl_matrix * LU, const gsl_permutation * p, gsl_matrix * inverse)
	/// call gsl_deallocate
	void matrix_inverse (gsl_matrix* inverse_result,int dimension_sqr);


	/// int gsl_matrix_transpose (gsl_matrix * m)
	// int gsl_matrix_transpose_memcpy (gsl_matrix * dest, const gsl_matrix * src)
	void matrix_transpose (gsl_matrix*  transposed);

	/// int gsl_matrix_add (gsl_matrix * a, const gsl_matrix * b)
	void matrix_add ( gsl_matrix * result, gsl_matrix* b);


	// Read Matrix from file: file user1 user2

	void LS_closed(CSR & Rat,Dense_matrix &Y,int K,float Lemda); 
	
	// generate Dense matrix from CSr formate
	void Generate_from_CSR(CSR & Rat);
	// generate CSR Matrix from Dense matrix
	void Generate_CSR(CSR & Rat,int r,int c,int v);
 
	// generate Dense matrix from file: file formate: node_id node_id 1 (which shows line from node_id to node_id)
	void Generate_from_File(std::ifstream &trfile,int r, int c);

	// Matrix multiplication for trust propogation for k times
	void Propogate_trust(Dense_matrix &T_rat,int k);
};

