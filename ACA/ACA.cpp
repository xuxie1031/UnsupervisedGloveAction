//============================================================================
// Name        : ACA.cpp
// Author      : XuXie
// Version     :
// Copyright   : XuXie@UCLA
// Description : C++, Ansi-style
//============================================================================

#include "ACA.h"

bool init_params(int &n, int &k, int &len_seg_pos, std::string data_name, int NC)
{
	struct stat info;
	if(stat("ACAbin", &info) != 0 || !S_ISDIR(info.st_mode))
		return false;

	std::stringstream ss;
	ss << NC << "_" << data_name;

	FILE* fid = fopen(("ACAbin/params_C_" + ss.str()).c_str(), "rb");
	if(!fid) return false;
	fread(&n, sizeof(int), 1, fid);
	fread(&k, sizeof(int), 1, fid);
	fread(&len_seg_pos, sizeof(int), 1, fid);
	fclose(fid);

	return true;
}

bool init_vars(int* Cids, int** seg_pos, double** Kb, int* mcs, double* thau_YYs, int n, int k, int len_seg_pos, std::string data_name, int NC)
{
	struct stat info;
	if(stat("ACAbin", &info) != 0 || !S_ISDIR(info.st_mode))
		return false;

	std::stringstream ss;
	ss << NC << "_" << data_name;

	FILE* fid = fopen(("ACAbin/thauYYs_C_" + ss.str()).c_str(),"rb");
	if(!fid) return false;
	fread(thau_YYs, sizeof(double), k, fid);
	fclose(fid);

	fid=fopen(("ACAbin/segpos_C_" + ss.str()).c_str(),"rb");
	if(!fid) return false;
	for(int i=0;i<len_seg_pos;i++)
		fread(seg_pos[i],sizeof(int),2,fid);
	fclose(fid);

	fid=fopen(("ACAbin/Cids_C_" + ss.str()).c_str(),"rb");
	if(!fid) return false;
	fread(Cids,sizeof(int),n,fid);
	fclose(fid);

	fid=fopen(("ACAbin/Kb_C_" + ss.str()).c_str(),"rb");
	if(!fid) return false;
	for(int i=0;i<n;i++)
		fread(Kb[i],sizeof(double),n,fid);
	fclose(fid);

	fid=fopen(("ACAbin/mcs_C_" + ss.str()).c_str(),"rb");
	if(!fid) return false;
	fread(mcs,sizeof(int),k,fid);
	fclose(fid);

	return true;
}

void minimum_dist(double &argdist, int &argc, int** seg_pos, double*** U, double** Kb, int* labels, int* mcs, double* thau_YYs, int v, int nv, int k, int len_seg_pos)
{
	double thau_XYs[k];
	memset(thau_XYs,0.0,sizeof(thau_XYs));
	argdist=1e5;
	argc=-1;

	for(int i=0;i<len_seg_pos;i++)
	{
		int idy=seg_pos[i][0];
		int len_segy=seg_pos[i][1]-seg_pos[i][0]+1;
		if(nv==0||v==0)
		{
			for(int j=idy+1;j<idy+len_segy;j++)
			{
				U[nv][v%nmax][j]=U[nv][v%nmax][j-1]+Kb[v][j];
			}
		}
		else
		{
			for(int j=idy+1;j<idy+len_segy;j++)
			{
				U[nv][v%nmax][j]=U[nv][v%nmax][j-1]+Kb[v][j];
				U[nv][v%nmax][j]=std::max(U[nv][v%nmax][j],U[nv-1][(v-1)%nmax][j-1]+2*Kb[v][j]);
				U[nv][v%nmax][j]=std::max(U[nv][v%nmax][j],U[nv-1][(v-1)%nmax][j]+Kb[v][j]);
			}
		}
		thau_XYs[labels[idy]]+=U[nv][v%nmax][idy+len_segy-1]/(nv+1+len_segy);
	}

	double sum_thau_XYs=0.0;
	for(int c=0;c<k;c++)
	{
		thau_XYs[c]=thau_XYs[c]*2/mcs[c];
		sum_thau_XYs+=thau_XYs[c];
	}
	for(int c=0;c<k;c++)
		thau_XYs[c]/=sum_thau_XYs;

	for(int c=0;c<k;c++)
	{
		double dist=sqrt(1.0-thau_XYs[c]+thau_YYs[c]);
		if(dist<argdist)
		{
			argdist=dist;
			argc=c;
		}
	}
}

void JACA(int** seg_pos, double** Kb, int* Cids, int* mcs, double* thau_YYs, int n, int k, int len_seg_pos)
{
	double*** U;
	Init3DArray<double>(U,nmax,nmax,n);
	double* J=new double[n];
	int* g=new int[n];
	int* iv=new int[n];

	printf("\nJdp now...need some time...\n");
	for(int v=0;v<n;v++)
	{
		J[v]=(double)INT_MAX;
		for(int i=0;i<len_seg_pos;i++)
		{
			int idy=seg_pos[i][0];
			U[0][v%nmax][idy]=2*Kb[v][idy];
		}

		for(int nv=0;nv<nmax;nv++)
		{
			int i=v-nv;
			if(i>=0)
			{
				double argdist;
				int argc;

				minimum_dist(argdist,argc,seg_pos,U,Kb,Cids,mcs,thau_YYs,v,nv,k,len_seg_pos);
				double tmpJ;
				if(i==0)
					tmpJ=argdist;
				else
					tmpJ=J[i-1]+argdist;
				if(tmpJ<J[v])
				{
					J[v]=tmpJ;
					g[v]=argc;
					iv[v]=i;
				}
			}
		}
	}

	int v=n-1;
	while(v>=0)
	{
		for(int i=iv[v];i<v+1;i++)
			Cids[i]=g[v];
		v=iv[v]-1;
	}

	Delete3DArray<double>(U,nmax,nmax);
	delete[] J;
	delete[] g;
	delete[] iv;
}

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		printf("arguments should be: \'number of clusters(integer)\' \'associated data name(string)\'\neg. ./ACA 9 sample_hand_data\n");
		return -1;
	}

	int NC;
	std::stringstream ss;
	ss << std::string(argv[1]);							// parse number of clusters
	ss >> NC;

	std::string data_name(argv[2]);						// parse associated data name


	// std::string data_name("sample_hand_data");			// associated data name
	// int NC = 9;											// number of cluster we're interested to apply ACA on
	int n, k, len_seg_pos;								// scalar parameters

	// start ACA optimization
	printf("initializing params for ACA...\n");
	if(init_params(n, k, len_seg_pos, data_name, NC))
	{
		printf("initializing vars for ACA...\n");

		// binaries initialization
		double** Kb;
		Init2DArray<double>(Kb,n,n);

		int* Cids=new int[n];

		int** seg_pos;
		Init2DArray<int>(seg_pos,len_seg_pos,2);

		int* mcs=new int[k];

		double* thau_YYs=new double[k];

		if(init_vars(Cids, seg_pos, Kb, mcs, thau_YYs, n, k, len_seg_pos, data_name, NC))
		{
			printf("ACA optimization...\n");

			JACA(seg_pos, Kb, Cids, mcs, thau_YYs, n, k, len_seg_pos);

			std::stringstream ss;
			ss << "ACAbin/ACA_Cids_C_" << NC << "_" << data_name;
			FILE* fid = fopen(ss.str().c_str(), "wb");
			fwrite(Cids, sizeof(int), n, fid);
			fclose(fid);
		}

		// binaries release
		Delete2DArray<double>(Kb,n);
		Delete2DArray<int>(seg_pos,len_seg_pos);
		delete[] Cids;
		delete[] mcs;
		delete[] thau_YYs;
	}

	return 0;
}
