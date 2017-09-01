#include <cmath>
#include "runtimecounter.h"
#include "Dense_Matrix.h"


using namespace std;

//int N_user = 49290;
//int N_item = 139738;
//int N_reviews = 664824;

vector<int> list,valid;
// validation + testing  = 30 % = 199447
// validation count 15 %
int N_val = 99724;
// testing count 15 %
int N_test = 99723;

int N_user = 6;
int N_item = 6;
int N_reviews = 19;

int K;
int N;
Runtimecounter rt,rt1;



float Get_objVal(CSR &Mat,Dense_matrix &Trust_Rat,Dense_matrix &P,Dense_matrix &Q,float Lemda1,float Lemda2, float alpha)
{
	float ret;
	float norm1 = P.get_frobenious_norm ();
	float norm2 = Q.get_frobenious_norm ();
	float normP = norm1*norm1;
	float normQ = norm2*norm2;
	float sum1 = 0.0;
	float sum2 = 0.0;

	for(int i=0;i < Mat.R ; i++)
	{
		int st = Mat.row_ptr[i];
		int end = Mat.row_ptr[i+1];
		int count = Mat.row_ptr[i+1] - st;
		gsl_matrix_view a = gsl_matrix_submatrix (P.data_,i,0,1,P.ncol_);
//		pi = P[i][:];
		int ind1 = 0;
		int counter = 0;
		for(int j= 0 ; j < Mat.C;j++)
		{
//			int col = Mat.col_ind[j];
			int col = j;
			gsl_matrix_view b = gsl_matrix_submatrix (Q.data_,col,0,1,Q.ncol_);
			
			// Transpose
			Dense_matrix d_temp (Q.ncol_,1);						// Q.ncol_ = K value
//			d_temp.initialize_matrix(0.0);
			gsl_matrix_set_zero(d_temp.data_);
			gsl_matrix_transpose_memcpy (d_temp.data_, &b.matrix);
			//multiplication result matrix
			Dense_matrix dx (1,1);
			gsl_matrix_set_zero(dx.data_);
//			dx.initialize_matrix(0.0);
			// multiplication (dot product)
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0,&a.matrix,d_temp.data_,0.0,dx.data_);


			if( ind1 < count && j == Mat.col_ind[st+ind1])
			{
				sum1 += (Mat.val[st+ind1] - gsl_matrix_get(dx.data_,0,0))* (Mat.val[st+ind1] - gsl_matrix_get(dx.data_,0,0));
				ind1++;
			}
			else
			{
//				sum2 += (gsl_matrix_get(Trust_Rat.data_,i,j) - gsl_matrix_get(dx.data_,0,0))* (gsl_matrix_get(Trust_Rat.data_,i,j) - gsl_matrix_get(dx.data_,0,0));
			}

//			qj = Q[b][:];
//			v = np.dot(pi,qj.T);
//			res = (Mat.val[j]-v)**2;
//			sum1 = sum1 + res;
		}
	}
	ret = sum1 + alpha*sum2  + Lemda1*normP + Lemda2*normQ;
	ret = ret/2.0;
	return ret;
}
double Get_MSE(CSR &Test,Dense_matrix &P,Dense_matrix &Q,float meanR,float Max,float Min)
{
	float MSE = 0.0;
	for(int i=0; i< Test.R; i++)
	{
		int st = Test.row_ptr[i];
		int end = Test.row_ptr[i+1];
		gsl_matrix_view a = gsl_matrix_submatrix (P.data_,i,0,1,P.ncol_);
		for(int j = st; j < end; j++)
		{
			int col = Test.col_ind[j];
			gsl_matrix_view b = gsl_matrix_submatrix (Q.data_,col,0,1,Q.ncol_);
			
			// Transpose
			Dense_matrix d_temp (Q.ncol_,1);						// Q.ncol_ = K value
			gsl_matrix_set_zero(d_temp.data_);
			gsl_matrix_transpose_memcpy (d_temp.data_, &b.matrix);

			//multiplication result matrix
			Dense_matrix dx (1,1);
			gsl_matrix_set_zero(dx.data_);

			// multiplication (dot product)
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0,&a.matrix,d_temp.data_,0.0,dx.data_);
//			Pred = np.dot(pi,qj.T);

			double Pred = gsl_matrix_get(dx.data_,0,0) + (double)meanR;
			if (Pred > Max)
				Pred = Max;
			if (Pred < Min)
				Pred = Min;
//			cout << i << " - " << j << " -> " << Pred << endl;
			double SE = ((double)Test.val[j]-Pred)*((double)Test.val[j]-Pred);
			MSE = MSE + SE;
		}
	}
	MSE = MSE / Test.V;

	return MSE;
}
int main(int argc, char **argv)
{
	char *trainfile, *testfile, *trustfile, *K_Lfile, * outfile;
	if(argc != 3)
		cout<< "Give  training file name, testing file, K-N value file name and output file name";
	else
	{		
		trainfile = argv[1];
		testfile = argv[2];
//		K_Lfile = argv[3];
//		outfile = argv[4];
	}


//---------------------------------- Part 1 --------------------------------------------------
	//training file preocessing
	ifstream trfile,tstfile;
	CSR Mat,Trans,Val,Test;

//	Generate testing samples
	//Generate_testing();

	cout << "file read & loaded: "<<endl;
	cout.flush();
	trfile.open (trainfile);
	if(trfile.is_open())
		Mat.loadTOCSR(trfile);
//		loadTOCSR(trfile);
	else
		cout << "cant open file";
	cout << "Done." <<endl;


	//sort can be done by this method
	Dense_matrix tmp;
	tmp.Generate_from_CSR(Mat);
	Mat.free();
	tmp.Generate_CSR(Mat,Mat.R,Mat.C,Mat.V);

//	Mat.disp();

	float meanR,Min, Max;
	Mat.Get_mean(meanR,Max,Min);
	Mat.Adjust_by_mean(meanR);


	cout << "matrix transposed: "<<endl;
	cout.flush();
	rt.start();
	transpose(Mat,Trans);
	cout << "Done."<<endl;


	int K,maxIter;
	float Lemda1,Lemda2, epsilon;
	float alpha, beta1, beta2;

	K = 20;
	maxIter = 100;
	Lemda1 = 3.0;
	Lemda2 = 3.0;
	epsilon = 0.00001;
	alpha = 0.0;

	Dense_matrix P(Mat.R,K);
	P.initialize_matrix (1.0/K);
	
	Dense_matrix Q(Mat.C,K);
	Q.initialize_matrix (1.0/K);

	cout << "Test file read & loaded: "<<endl;
	tstfile.open (testfile);
	if(tstfile.is_open())
		Test.loadTOCSR(tstfile);
//		loadTOCSR(trfile);
	else
		cout << "cant open file";
	cout << "Done." <<endl;

	int t =0;
	float cur;
	float pre_obj = Get_objVal(Mat,Trust_Rating,P,Q,Lemda1,Lemda2,alpha);

	cout << pre_obj;

	double min_MSE = 99999.0;
	while(t < maxIter)
	{
		t = t +1;
		P.LS_closed(Mat,Q,K,Lemda1);
		Q.LS_closed(Trans,P,K,Lemda1);

//		correct LS_closed function
/*		for(int i =0 ;i < 10 ; i++)
		{
			cout << i << " -> " << gsl_matrix_get(Q.data_,i,0);
			cout << endl;
		}
//*/
		cur = Get_objVal(Mat,Trust_Rating,P,Q,Lemda1,Lemda2,alpha);

		float check = abs(pre_obj-cur);
		check = check/pre_obj;
		if (check < epsilon)
			break;
		double M = Get_MSE(Test,P,Q,meanR,Max,Min);
		if (M < min_MSE)
			min_MSE = M;
		printf("%d -> %f  -  %f\n",t,cur,M);
	}
//--------------------------------- Part 2 ----------------------------------------------------

//~~~~~~~~~~~~~  Part 2(a)

	printf("Min MSE = %f\n",min_MSE);


	Mat.free();
//	Sim_mat.free();
	Test.free();
	return 0;
}
