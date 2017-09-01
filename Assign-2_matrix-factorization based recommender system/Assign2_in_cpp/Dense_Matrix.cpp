#include "Dense_Matrix.h"

using namespace std;

/// Constructor
Dense_matrix::Dense_matrix (){}
Dense_matrix::Dense_matrix (size_t nrow, size_t ncol)
: nrow_(nrow),
ncol_(ncol),
data_(gsl_matrix_alloc(nrow,ncol))
{
}
void Dense_matrix::init_alloc (size_t nrow, size_t ncol)
{
	nrow_ = nrow;
	ncol_ = ncol,
	data_ = gsl_matrix_alloc(nrow,ncol);
}

/// Destructor
Dense_matrix::~Dense_matrix ()
{
	gsl_matrix_free (data_);
}
/// Initialize the matrix with a particular value
void Dense_matrix::initialize_matrix (double val)
{
	gsl_matrix_set_all (data_,val);
}

// 2nd order forbious norm
float Dense_matrix::get_frobenious_norm ()
{
	float frob =0.0;
	float val;
	for (int i=0; i<nrow_; i++)
	{
		for (int j=0; j<ncol_; j++)
		{
			val = gsl_matrix_get(data_,i,j);
			frob = frob+ (val*val);
		}   
	}
	return sqrt(frob);
}

// int gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C)
void Dense_matrix::matrix_multiply (gsl_matrix* b, gsl_matrix* result)
{
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, this->data_, b,0.0,result);
}	

/// gsl_matrix invert
/// gsl_permutation * gsl_permutation_alloc (size_t n)
/// void gsl_permutation_free (gsl_permutation * p)
/// int gsl_linalg_LU_decomp (gsl_matrix * A, gsl_permutation * p, int * signum)
/// int gsl_linalg_LU_invert (const gsl_matrix * LU, const gsl_permutation * p, gsl_matrix * inverse)
/// call gsl_deallocate

void Dense_matrix::matrix_inverse (gsl_matrix* inverse_result,int dimension_sqr)
{
	int s;
	Dense_matrix d1(this->nrow_,this->ncol_);
	gsl_matrix_memcpy (d1.data_, this->data_);
	gsl_permutation * perm = gsl_permutation_alloc (dimension_sqr);
	gsl_linalg_LU_decomp (d1.data_, perm, &s);
	gsl_linalg_LU_invert (d1.data_, perm, inverse_result);
	gsl_permutation_free (perm);
}

/// int gsl_matrix_transpose (gsl_matrix * m)
// int gsl_matrix_transpose_memcpy (gsl_matrix * dest, const gsl_matrix * src)
void Dense_matrix::matrix_transpose (gsl_matrix*  transposed)
{
	gsl_matrix_transpose_memcpy (transposed,this->data_);
}
/// int gsl_matrix_add (gsl_matrix * a, const gsl_matrix * b)
void Dense_matrix::matrix_add ( gsl_matrix * result, gsl_matrix* b)
{
	gsl_matrix_memcpy (result, this->data_);
	gsl_matrix_add (result, b);
}

void Dense_matrix::LS_closed(CSR &Rat,Dense_matrix &Y,int K,float Lemda)
{

	Dense_matrix* outer_list = new Dense_matrix[Rat.C];

	Dense_matrix sum1(1,K),result_temp(1,K), inverse_result(K,K);


//		d1.initialize_matrix(0.0); 
//		result_temp.initialize_matrix(0.0);
//		inverse_result.initialize_matrix(0.0);

	Dense_matrix matrix(K,K);
//		matrix.initialize_matrix(0.0);


	Dense_matrix I(K,K);
	gsl_matrix_set_identity(I.data_);
	// multiply Lemda with Identity
	for(int i =0; i < K ; i++)
		for(int j =0 ; j < K ; j++)
			if(i==j)
				gsl_matrix_set(I.data_,i,j,gsl_matrix_get(I.data_,i,j) * Lemda);

	Dense_matrix d_temp (K,1);
	Dense_matrix myq(1,K); 

	for(int i =0 ; i < Y.nrow_ ; i++)
	{

		gsl_matrix_set_zero(myq.data_);
//			myq.initialize_matrix(0.0);
		for(int j =0 ; j < Y.ncol_ ; j++)					// ncol_= K
		{
			double val = gsl_matrix_get(Y.data_,i,j);
			gsl_matrix_set(myq.data_,0,j,val);
		}

//			d_temp.initialize_matrix(0.0);
		gsl_matrix_set_zero(d_temp.data_);
		gsl_matrix_transpose_memcpy (d_temp.data_,myq.data_);

		outer_list[i].init_alloc(K,K);

		// do the multiplication
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, d_temp.data_,myq.data_,0.0,outer_list[i].data_);

//			t = Y[i][:];
//			mul = np.outer(t,t);
//		outer_list.push_back(matrix_temp);
	}
	cout << "Here"<<endl;
//*/
	for(int i=0; i < Rat.R; i++)
	{
		int st = Rat.row_ptr[i];
		int end = Rat.row_ptr[i+1];
		int count = end-st;
		gsl_matrix_set_zero(sum1.data_);
		gsl_matrix_set_zero(matrix.data_);
//			sum1 = np.zeros((1,K), dtype=np.float);
//			sum2 = np.zeros((K,K), dtype=np.float);

		int ind1 = 0;
		for(int j = 0; j <Rat.C; j++)
		{
//			int col = Rat.col_ind[j];
			int col = j;
			float value = Rat.val[st+ind1];

			Dense_matrix myq(1,K); 
			gsl_matrix_set_zero(myq.data_);
			Dense_matrix myq_R(1,K); 
			gsl_matrix_set_zero(myq.data_);
//			myq.initialize_matrix(0.0);
			for(int jj =0 ; jj < Y.ncol_ ; jj++)					// ncol_= K
			{
				double val = gsl_matrix_get(Y.data_,col,jj);
				gsl_matrix_set(myq_R.data_,0,jj,(double)value * val);
				gsl_matrix_set(myq.data_,0,jj,val);
			}

		
//	#			s2 = np.outer(t,t);
/*
			gsl_matrix_set_zero(d_temp.data_);
			gsl_matrix_transpose_memcpy (d_temp.data_,myq.data_);

//			matrix_temp.initialize_matrix(0.0); 
			Dense_matrix matrix_temp(K,K);
			gsl_matrix_set_zero(matrix_temp.data_);

			// do the multiplication
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, d_temp.data_,myq.data_,0.0,matrix_temp.data_);
*/
//				sum2 = sum2 + s2;



			if( ind1 < count && j == Rat.col_ind[st+ind1])
			{
				gsl_matrix_add (sum1.data_, myq_R.data_);
//				s1 = Y[b][:] * X.val[j];
//				sum1 = sum1 + s1;

				gsl_matrix_add (matrix.data_, outer_list[col].data_);
//				matrix = matrix + matrix_temp;

				ind1++;
			}
			else
			{
			}


		}
		
		gsl_matrix_add (matrix.data_, I.data_);
//			r = I + sum2;

		matrix.matrix_inverse(inverse_result.data_,K);
//			secondT = np.linalg.inv(r);

		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, sum1.data_,inverse_result.data_,0.0,result_temp.data_);
//			res = np.dot(sum1,secondT);

		// set new row values
		for(int j =0;j < K; j++)
		{
			gsl_matrix_set(this->data_,i,j,gsl_matrix_get(result_temp.data_,0,j));
//			if(i < 5 && j < 20)
//				cout << gsl_matrix_get(result_temp.data_,0,j) << " - ";
		}


		gsl_matrix_set_zero(result_temp.data_);
		gsl_matrix_set_zero(inverse_result.data_);

	}

}

// generate Dense matrix from CSR formate
void Dense_matrix::Generate_from_CSR(CSR & Rat)
{
	this->init_alloc(Rat.R,Rat.C);
	gsl_matrix_set_zero(this->data_);
	for(int i =0 ;i < Rat.R; i++)
	{
		for(int j = Rat.row_ptr[i]; j < Rat.row_ptr[i+1]; j++)
		{
			gsl_matrix_set(this->data_,i,Rat.col_ind[j],Rat.val[j]);
		}
	}
}

// generate CSR Matrix from Dense matrix
void Dense_matrix::Generate_CSR(CSR & Rat,int r,int c,int v)
{
	
	Rat.init(r,c,v);
	int counter = 0;
	Rat.row_ptr[0] = 0;

	for(int i =0 ;i < r; i++)
	{
		for(int j = 0; j < c; j++)
		{
			float t = gsl_matrix_get(this->data_,i,j);
			if(t != 0.0)
			{
				Rat.col_ind[counter] = j;
				Rat.val[counter] = t;
				counter++; 
			}
		}
		Rat.row_ptr[i+1] = counter;
	}

//	Rat.disp();
}
// generate Dense matrix from file: file formate: node_id node_id 1 (which shows line from node_id to node_id)
void Dense_matrix::Generate_from_File(std::ifstream &trfile,int r, int c)
{
	this->init_alloc(r,c);
	gsl_matrix_set_zero(this->data_);


	string s;				//read line
	string delimiter = " ";


	while(true)
	{
		getline(trfile,s);
		if(s.size() == 0)
			break;
		string token;
		int i =0;
		int r,c,v;
		stringstream os(s);
		while(os >> token)
		{
			istringstream buffer(token);
			size_t t;
			buffer >> t;  
			if(i ==0)
				r = t;
			else if(i == 1)
				c = t;
			else
				v = t;
	//			printf("%d\n",t);
			i++;
		}
		gsl_matrix_set(this->data_,r-1,c-1,v);
	}

}

// Matrix multiplication for trust propogation for k times
void Dense_matrix::Propogate_trust(Dense_matrix &T_rat,int k)
{
	for(int i =0; i < k ; i++)
	{
		gsl_matrix * result = gsl_matrix_alloc(this->nrow_,this->ncol_);
		// do the multiplication
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, this->data_,T_rat.data_,0.0,result);

		// reset current matrix values to resultant values
		gsl_matrix_set_zero(this->data_);
		gsl_matrix_memcpy (this->data_,result);
	}
}
