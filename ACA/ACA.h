//============================================================================
// Name        : ACA.h
// Author      : XuXie
// Version     :
// Copyright   : XuXie@UCLA
// Description : C++, Ansi-style
//============================================================================


#include <stdio.h>
#include <string.h>
#include <sstream>
#include <math.h>
#include <limits.h>
#include <sys/stat.h>

#define nmax 100

template <typename T>
void Init2DArray(T** &arr, int d1, int d2)
{
	arr=new T*[d1];
	for(int i=0;i<d1;i++)
	{
		arr[i]=new T[d2];
	}
}

template <typename T>
void Delete2DArray(T** &arr, int d1)
{
	for(int i=0;i<d1;i++)
	{
		delete[] arr[i];
	}
}

template <typename T>
void Init3DArray(T*** &arr, int d1, int d2, int d3)
{
	arr=new T**[d1];
	for(int i=0;i<d1;i++)
	{
		arr[i]=new T*[d2];
		for(int j=0;j<d2;j++)
		{
			arr[i][j]=new T[d3];
		}
	}
}

template <typename T>
void Delete3DArray(T*** &arr, int d1, int d2)
{
	for(int i=0;i<d1;i++)
	{
		for(int j=0;j<d2;j++)
		{
			delete[] arr[i][j];
		}
	}
}


// load necessary parameters
bool init_params(int &n, int &k, int &len_seg_pos, std::string data_name, int NC);

// load variables to initialize ACA procedure
bool init_vars(int* Cids, int** seg_pos, double** Kb, int* mcs, double* thau_YYs, int n, int k, int len_seg_pos, std::string data_name, int NC);

// calculate the minimum dist
void minimum_dist(double &argdist, int &argc, int** seg_pos, double*** U, double** Kb, int* labels, int* mcs, double* thau_YYs, int v, int nv, int k, int len_seg_pos);

// ACA main entry
void JACA(int** seg_pos, double** Kb, int* Cids, int* mcs, double* thau_YYs, int n, int k, int len_seg_pos);