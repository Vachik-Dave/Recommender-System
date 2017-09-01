#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <vector>
#include "runtimecounter.h"

using namespace std;

int K;
int N;
Runtimecounter rt,rt1;

class CSR
{
public:
	int R,C,V;			// R - row size, C - column size, V - Nonzero size
	int *col_ind, *row_ptr;
	float * val;

	CSR(){};
	CSR(int r,int c,int v)
	{
		R = r; C = c; V = v;
	}
	void init(int r,int c,int v)
	{
		R = r; C = c; V = v;

		val = new float[V];
		col_ind = new int[V];
		row_ptr = new int[R+1];
	}
	float get_val(int r,int c)
	{
		for(int i =row_ptr[r] ; i < row_ptr[r+1] ; i++)
		{
			if(col_ind[i] == c)
				return val[i];
			if(col_ind[i] > c)
				return (0.0);
		}
		return (0.0);		
	}
	void disp()
	{
		int i;
		for(i =0;i < V;i++)
		{
			printf("%3d - ",col_ind[i]);
		}
		cout << endl;
		for(i =0;i < V;i++)
		{
			printf("%3f - ",val[i]);
		}
		cout <<endl;
		
		for(i=0;i<R+1;i++)
			printf("%3d - ",row_ptr[i]);
		cout << endl;
	}
	// Read given file and convert data into CSR format
	void loadTOCSR(ifstream &trfile)
	{
		string s;				//read line
		string delimiter = " ";
	
		getline(trfile,s);

		string token;
		int i =0;
		stringstream os(s);
		while(os >> token)
		{
			istringstream buffer(token);
			int t;
			buffer >> t;  
			if(i ==0)
				R = t;
			else if(i == 1)
				C = t;
			else
				V = t;
//			printf("%d\n",t);
			i++;

		}

		init(R,C,V);

		int count =0;
		for(i=0; i < R; i++)
		{
			row_ptr[i] = count;

			s.clear();
			getline(trfile,s);
			int j =0;
			stringstream os(s);
			while(os >> token)
			{
				istringstream buffer(token);
				if(j%2 ==0)
				{
					int t;
					buffer >> t;  
					if(t > C)
						cout << "Error: get column number higher than C(size)"<< endl;
					col_ind[count] = t-1;				// -1 because indexes in input file starts with 1 
				}
				else
				{
					float t;
					buffer >> t;  
					if(t==0)
						cout << "Error: get 0 input" <<endl;
					val[count] = t;
					count++;
				}
		//		printf("%d\n",t);
				j++;
			}	
			if(j%2 != 0)
				cout << "Error: tokening problem"<< endl;


		}
		if(count == V)
			row_ptr[R] = count;
		else
			cout << "Error: non zero count different" << endl;

//		disp();
	}

	void store_CSR(FILE* &store)
	{
//		FILE* store = fopen(filename,"w+");
		fprintf(store,"%d %d %d",R,C,V);
		int r = 0;
		for(int i =0; i < V; i++)
		{
			if(i == row_ptr[r])
			{
				fprintf(store,"\n");
				r++;
			}
			fprintf(store,"%d %3.6f ",col_ind[i]+1,val[i]);
		}
	}
	void free()
	{
		delete[] val;
		delete[] col_ind;
		delete[] row_ptr;
	}

}Mat,Trans,Sim_mat,Test;
vector<pair<int,float> > *item_sim;
vector<int> * sim_sort;
vector<pair<int,float> > * Reco;


float get_cos_sim(int i,int j,CSR & M)
{

	int nrowi = M.row_ptr[i+1] - M.row_ptr[i];
	int nrowj = M.row_ptr[j+1] - M.row_ptr[j];

	int ni = 0, nj = 0;
	float cosine = 0, lengthi = 0, lengthj = 0;

	while (ni < nrowi && nj < nrowj)
	{

		int ci = M.row_ptr[i] + ni;
		int cj = M.row_ptr[j] + nj;

		if (M.col_ind[ci] == M.col_ind[cj])
		{
/*
			cosine  += M.val[ci] * M.val[cj];
			lengthi += M.val[ci] * M.val[ci];
			lengthj += M.val[cj] * M.val[cj];
*/
			cosine  += 1 * 1;			// binary values
			lengthi += 1 * 1;
			lengthj += 1 * 1;
		
			ni++;nj++;
		}
	    	else if (M.col_ind[ci] > M.col_ind[cj])
		{
//			lengthj += M.val[cj] * M.val[cj];

			lengthj += 1 * 1;
			nj ++;
	    	}
		else
		{
//			lengthi += M.val[ci] * M.val[ci];

			lengthi += 1 * 1;
			ni ++;
	    	}
	}
	while(ni < nrowi)
	{
//		lengthi += M.val[ci] * M.val[ci];

		lengthi += 1 * 1;
		ni ++;
	}
	while(nj < nrowj)
	{
//		lengthj += M.val[cj] * M.val[cj];

		lengthj += 1 * 1;
		nj ++;
	}
	
//	cout << i << " - " << j << " : " << cosine << " / " << lengthi * lengthj << endl;

	if (lengthi * lengthj != 0)
	    cosine /= sqrt(lengthi * lengthj);
	else
	    cosine = 0;

	return cosine;
}

void transpose()
{

	Trans.init(Mat.C,Mat.R,Mat.V);				//Rt = C; 	Ct = R;		Vt = V;

	int * row_counts2 = new int[Trans.R]();		//set all to zero 

	for(int i  =0; i < Trans.V; i++)
	{
		int t = Mat.col_ind[i];
		row_counts2[t]++;
	}

	Trans.row_ptr[0] =0;
	for(int i =1 ;i < Trans.R+1; i++)
	{
		Trans.row_ptr[i] = Trans.row_ptr[i-1] + row_counts2[i-1] ;
	}

	memset(row_counts2, 0, sizeof(int) * Trans.R);			// reinitialize counters to zero


	for (int i = 0; i < Mat.R; i ++)
	{
		for (int j = Mat.row_ptr[i]; j < Mat.row_ptr[i+1]; j ++)
		{
			int c = Mat.col_ind[j]; 

			if(row_counts2[c] > (Trans.row_ptr[c+1]-Trans.row_ptr[c]))
				cout << "Error: count proble in trnspose" << endl;

			Trans.col_ind[Trans.row_ptr[c] + row_counts2[c]] = i; 
			Trans.val[Trans.row_ptr[c] + row_counts2[c]] = Mat.val[j];  
			row_counts2[c] ++;  
	    	}
	}

	delete[] row_counts2;

//	Trans.disp();
}


inline bool myfunction1 (pair<int,float> i,pair<int,float> j) 
{ 
	return (i.second > j.second); 
}

void Build_CSR_mat(int count)
{
	Sim_mat.init(Trans.R,Trans.R,count);

	int c =0;
	int counter =0;
	Sim_mat.row_ptr[0] = 0;
	for(int i =0 ; i < Trans.R; i++)
	{
		counter =0;
		for(int j = 0; j < item_sim[i].size(); j++)
		{
			if(item_sim[i][j].first <= i)
				continue;

			Sim_mat.col_ind[c] = item_sim[i][j].first;
			Sim_mat.val[c] = item_sim[i][j].second;
			c++;
			counter++;
		}
		Sim_mat.row_ptr[i+1] = Sim_mat.row_ptr[i]+counter;
	}
	if(c != count)
		cout << "Error: in counting" <<endl;

//	Sim_mat.disp();
	
}

void write_to_file(FILE *ofile)
{
	fprintf(ofile,"K value = %d & N value = %d\n",K,N);
//	printf("K value = %d & N value = %d\n",K,N);
	for(int i =0 ; i < Mat.R; i++)
	{
		int siz = Reco[i].size();
		int limit = N;
		if(siz < N)
			limit = siz;
		fprintf(ofile,"%d:",i+1);								//user_id = i + 1 bcoz i starts with 0
//		printf("%d:",i+1);
		for(int j =0; j < limit; j++)
		{
			fprintf(ofile," %d %5.2f",Reco[i][j].first+1,Reco[i][j].second);			// item_id = index+1
//			printf("%d %5.2f",Reco[i][j].first+1,Reco[i][j].second);
		}

		fprintf(ofile,"\n");
//		printf("\n");
	}

}


int main(int argc, char **argv)
{
	char *trainfile, *testfile, *K_Nfile, * outfile;
	if(argc != 5)
		cout<< "Give  training file name, K-N value file, testing file name and output file name";
	else
	{		
		trainfile = argv[1];
		testfile = argv[2];
		K_Nfile = argv[3];
		outfile = argv[4];
	}


//---------------------------------- Part 1 --------------------------------------------------
	//training file preocessing
	ifstream trfile;
	cout << "file read & loaded: ";
	trfile.open (trainfile);
	if(trfile.is_open())
		Mat.loadTOCSR(trfile);
//		loadTOCSR(trfile);
	else
		cout << "cant open file";
	cout << "Done." <<endl;
	
	//K & N reading
	ifstream KNfile;
	KNfile.open (K_Nfile);

	string s;				//read line
	string delimiter = " ";

	vector<int> Klist;
	vector<int> Nlist;
	string token;
	

	// first row list of K values
	getline(KNfile,s);
	stringstream os(s);
	while(os >> token)
	{
		istringstream buffer(token);
		int t;
		buffer >> t;
		Klist.push_back(t);
	}
	// second row list of N values
	s.clear();
	getline(KNfile,s);
	stringstream os1(s);
	while(os1 >> token)
	{
		istringstream buffer(token);
		int t;
		buffer >> t;
		Nlist.push_back(t);
	}
	sort(Klist.begin(),Klist.end());
	sort(Nlist.begin(),Nlist.end());


//--------------------------------- Part 2(a) ----------------------------------------------------

//~~~~~~~~~~~~~  Part 2(a) <i>

	cout << "matrix transposed: ";
	rt.start();
	transpose();
	cout << "Done."<<endl;

//	item_sim = new vector<pair<int,float> > [Trans.R];

	float * * Simi_mat = new float* [Mat.C];
	for(int i =0 ; i < Mat.C; i++)
	{
		Simi_mat[i] = new float[Mat.C];
	}
	int count =0;


	cout << "similarities found: ";
//	cout << "here1";
	for(int i =0; i < Trans.R; i++)
	{
//	cout << i << " - ";
		for(int j = i+1; j < Trans.R; j++)				// j = i+1 b'coz cos_sim(i,j) = cos_sim(j,i)
		{
			float ret = get_cos_sim(i,j,Trans);
//			cout << i << " - " << j << " : " << ret << endl;
			if(ret > 0)
			{
//				item_sim[i].push_back(make_pair(j,ret));
//				item_sim[j].push_back(make_pair(i,ret));	
				Simi_mat[i][j] = ret;
				Simi_mat[j][i] = ret;
				count++;
			}
		}
	}
	cout << "Done."<<endl;
	rt.stop();						// convert & storing into CSR format should not be added in time
	float part_a1 = rt.GetRuntime();


	Trans.free();

//~~~~~~~~~~~~~  Part 2(a) <iv>	// sparce matrx representation of similarity
//	Build_CSR_mat(count);
/*	float * * Simi_mat = new float* [Mat.C];
	for(int i =0 ; i < Mat.C; i++)
	{
		Simi_mat[i] = new float[Mat.C];

		for(int j =0; j < item_sim[i].size();j++)
		{
			Simi_mat[i][item_sim[i][j].first] = item_sim[i][j].second;
		}
	}
*/
//~~~~~~~~~~~~~  Part 2(a) <iii>
	cout << "Sorted Similarities: ";		
	vector< pair<int,float> > item_sim;
	sim_sort = new vector<int> [Mat.C];						// stores similarity based sorted item_id for each item
	rt.start();
	for(int i = 0; i < Mat.C; i++)
	{
		item_sim.clear();
		for(int j =0; j < Mat.C; j++)
		{
			if(j != i || Simi_mat[i][j] != 0)
				item_sim.push_back(make_pair(j,Simi_mat[i][j]));
		}

		sort(item_sim.begin(),item_sim.end(),myfunction1);

		for(int j =0; j < item_sim.size(); j++)
		{
			sim_sort[i].push_back(item_sim[j].first);
		}
	}
	cout << "Done."<<endl;
//~~~~~~~~~~~~~  Part 2(a) <ii>
	rt.stop();
	float part_a2 = rt.GetRuntime();
	cout << "Part (a) timing = "<<part_a1+part_a2 <<" ms"<<endl;
//	FILE *file = fopen("Similarity_mat.txt","w");
//	Sim_mat.store_CSR(file);


//	delete[] item_sim;

//--------------------------------- Part 2(b) ----------------------------------------------------

	FILE *o_file = fopen(outfile, "w+");
	ifstream tstfile;
	tstfile.open (testfile);
	Test.loadTOCSR(tstfile);
	FILE* resfile = fopen("results.txt","w+");


	for(int T = 0; T < Klist.size(); T++)
	{
		K = Klist[T];
		cout << "K = " << K<<endl;

		// Memory assignment
		vector<int> * candidate = new vector<int>[Mat.R];			// 1 candidate list for each user

		// for Part 2(b) <ii> 
		Reco = new vector<pair<int,float> > [Mat.R];				//for each user list of recommendation



	//~~~~~~~~~~~~~  Part 2(b) <i>

		rt.start();
		cout << "Candidates found: ";
		for(int i =0; i < Mat.R; i++)
		{
	//		int count = Mat.row_ptr[i] - Mat.row_ptr[i+1];
			vector<int> v;
			for(int j = Mat.row_ptr[i] ; j < Mat.row_ptr[i+1]; j++)
				v.push_back(Mat.col_ind[j]);					// already rated list	

			for(int j = Mat.row_ptr[i] ; j < Mat.row_ptr[i+1]; j++)
			{
				int item_id = Mat.col_ind[j];
				int c =0;
				for(int k =0 ; k < sim_sort[item_id].size() ; k++)
				{
					int curr = sim_sort[item_id][k];
					c++;
					if(c >= K)
						break;
					if(find(v.begin(), v.end(), curr) != v.end())		// if next similar item is already rated
						continue;
					if(c < K)						// upto K most similar items
					{
						candidate[i].push_back(curr);
					}

				}
			}
			// remove duplicate
			vector<int>::iterator it;
			sort(candidate[i].begin(),candidate[i].end());
			it = unique(candidate[i].begin(),candidate[i].end());
			candidate[i].erase(it, candidate[i].end());
		}
		cout << "Done." <<endl;
	//~~~~~~~~~~~~~  Part 2(b) <ii>

		cout << "Recommendation list: ";
		float a=0.0,b=0.0,c=0.0;
		for(int i =0; i < Mat.R; i++)
		{
			vector<int> v;
			rt1.start();
			for(int j = Mat.row_ptr[i] ; j < Mat.row_ptr[i+1]; j++)
				v.push_back(Mat.col_ind[j]);
			rt1.stop();
			a = a+rt1.GetRuntime();

			rt1.start();
			int counter =0;
			for(int j =0 ; j < candidate[i].size(); j++)
			{
				int curr = candidate[i][j];
//				float s = get_score(curr,v);
				float sum = 0.0;
				for(int temp= 0;temp<v.size();temp++)
				{	counter++;
					int curr1 = v[temp];
//too slow				float tt = get_cos_sim(curr,curr1,Trans);
//					sum = sum + tt;

					float tt = Simi_mat[curr][curr1];
					sum = sum + tt;


/*	very slow
					if(curr1 > curr)
					{
						float tt = Sim_mat.get_val(curr,curr1);
						if(tt > 0)
							sum = sum + tt;
					}
					else if(curr1 < curr)
					{
						float tt = Sim_mat.get_val(curr1,curr);
						if(tt > 0)
							sum = sum + tt;
					}
					else
						cout << "Error: get_score() same values shouldnt be there" <<endl;
//*/
				}
				float s = sum;
				Reco[i].push_back(make_pair(curr,s));
			}
			if(counter != candidate[i].size() * v.size())
				cout << "counter = " << counter << " multiplication = "<< candidate[i].size() * v.size() <<endl;
			rt1.stop();
			b = b+rt1.GetRuntime();

			rt1.start();	
			sort(Reco[i].begin(),Reco[i].end(),myfunction1);
			rt1.stop();
			c = c+rt1.GetRuntime();

		}
		cout<< "Done."<<endl;
		cout << " a = " << a<< " b = " << b<< " c  = " << c <<endl;
		rt.stop();									//file writing time should not be consider
		float part_b1 = rt.GetRuntime();
		cout <<"For K =  " <<K<< " \t Part (b) timing = "<< part_b1<<" ms"<<endl;					
		
	//~~~~~~~~~~~~~  Part 2(b) <iii>
		for(int TT = 0; TT < Nlist.size();TT++)
		{
			N = Nlist[TT];

			write_to_file(o_file);
		}

	//--------------------------------- Part 2(c) ----------------------------------------------------


		for(int TT = 0; TT < Nlist.size();TT++)
		{
			rt.start();
			N = Nlist[TT];

			float HR=0.0,ARHR=0.0;
			for(int i =0; i < Mat.R; i++)
			{
				for(int j =0; j < N; j++)
				{
					if(j >= Reco[i].size())
						break;
					if(Reco[i][j].first == Test.col_ind[Test.row_ptr[i]])
					{

						HR++;
						ARHR = ARHR + (1.0/(float)(j+1));		// j could be 0 & rank starts with 1 so j+1
						break;
					}
				}
			}

			HR = HR / Mat.R;
			ARHR = ARHR / Mat.R; 

			fprintf(resfile,"%d %d %f %f\n",K,N,HR,ARHR);
			printf("%d %d %f %f",K,N,HR,ARHR);

			rt.stop();
			float part_b2 = rt.GetRuntime();
			//cout << " K = "<< K <<" & N = "<<N <<endl;
			cout << "\t Part (c) timing = "<< part_b2<<" ms"<<endl;

		}//for loop for N list

		delete[] candidate;
		delete[] Reco;
	}// for loop for K list


	Mat.free();
//	Sim_mat.free();
	Test.free();
	return 0;
}
