#include <fstream>
#include <sstream>
#include <string.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h> 
#include <algorithm>
#include <math.h>
#include <vector>
#include <iterator>


class CSR
{
public:
	int R,C,V;			// R - row size, C - column size, V - Nonzero size
	int *col_ind, *row_ptr;
	float * val;
//	int * val;

	// Member functions
//constructors
	CSR();
	CSR(int r,int c,int v);
//distructor
	~CSR();
// initializer to allocate memory
	void init(int r,int c,float v);
// Accessor
	float get_val(int r,int c);
// Data Statistics
	void Get_mean(float & mean,float & Max,float &Min);
// Manipulation for normalization
	void Adjust_by_mean(float mean);
// Display
	void disp();
// Read given file and convert data into CSR format
	void loadTOCSR(std::ifstream &trfile);
// Store data into a file in specific formate described in Recommended System 1st Assignment
	void store_CSR(FILE* &);
// Freeing memory intentionally
	void free();
};

