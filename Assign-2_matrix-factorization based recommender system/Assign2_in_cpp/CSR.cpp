#include "CSR.h"

using namespace std;

//constructors
CSR::CSR(){};
CSR::CSR(int r,int c,int v)
{
	R = r; C = c; V = v;
}
//distructor
CSR::~CSR()
{
	delete[] val;
	delete[] col_ind;
	delete[] row_ptr;
}
// initializer to allocate memory
void CSR::init(int r,int c,float v)
{
	R = r; C = c; V = v;

	val = new float[V];
	col_ind = new int[V];
	row_ptr = new int[R+1];
}
// Accessor
float CSR::get_val(int r,int c)
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
// Data Statistics
void CSR::Get_mean(float & mean,float & Max,float &Min)
{
	float sum1 = 0.0;
	Max = 0.0;
	Min = 99999.0;
	for(int i = 0; i < this->V; i++)
	{
		sum1 = sum1 + val[i];
		if(Max < val[i])
			Max = val[i];
		if(Min > val[i])
			Min = val[i];
	}
	mean = sum1 / this->V;
}
// Manipulation for normalization
void CSR::Adjust_by_mean(float mean)
{
	for (int i =0 ; i < V ;i++)
		val[i] = val[i] - mean;
}
//Display
void CSR::disp()
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
void CSR::loadTOCSR(std::ifstream &trfile)
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
// Store data into a file in specific formate described in Recommended System 1st Assignment
void CSR::store_CSR(FILE* &store)
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
		fprintf(store,"%d %f ",col_ind[i]+1,val[i]);
	}
}
// Freeing memory intentionally
void CSR::free()
{
	delete[] val;
	delete[] col_ind;
	delete[] row_ptr;
}

