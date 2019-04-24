#pragma once

#include "Headers.h"

void create_Bmat_bimodal(std::vector <float>& Bmat);

void update_Bmat_bimodal(std::vector <float>& Bmat, double del);

void create_Bmat_gaussian(std::vector <float>& Bmat, double del);

void create_Exmat(std::vector < std::vector <int> >& Exmat);

void create_CapacityMat(std::vector < std::vector <float> >& CapacityMat, std::vector < std::vector <int> >& Exmat);

void create_Wmat(std::vector < float >& Wmat, std::vector < std::vector <float> >& CapacityMat, std::vector <float >& BMat, double del);

void create_Augumented_CapacityMat(std::vector < float >& Wmat, std::vector < std::vector <float> >& CapacityMat);

void create_Residual_graph(std::vector < std::vector <float> >& CapacityMat, std::vector < std::vector <float> >& flow);

void printMatrix(std::vector <std::vector <int> >  & M, int len);

void printMatrix(std::vector <std::vector <float> >  & M, int len);

void printMatrix(std::vector<float> & M, int len);

void printMatrix(std::vector<int> & M, int len);

void savedata(std::vector< int > & Mat, int l, double del, std::string s, double omega, int rep);

void savedata(std::vector< std::vector< float > > & Mat, int l, double del, std::string s);

void savedata(std::vector< float > & Mat, int l, double del, std::string s, double omega, int rep);

void dfs(std::vector <std::vector <float> >& F, int s, std::vector <int>& visited);

void getSi(std::vector <int>& visited, int len);

double Hamiltonian(double maxflow, std::vector<std::vector<float>>& CapacityMat);

double Magnetization(std::vector <int>& visited);

double get_q(std::vector<int>& occ_init, std::vector<int>& occ);

void transform_phi(std::vector<float> &phi_in, std::vector<std::vector<float>> &phi_out);

void transform_si(std::vector<int> &si_in, std::vector<std::vector<int>> &si_out);

double Hamiltonian_formula(std::vector<float> &phi, std::vector<int> &si);

// this func re- intializes all the matrices with zero.
void Reinit(std::vector < std::vector <int> >&   sqlat0,
	std::vector < std::vector <int> >& sqlat1,
	std::vector < std::vector <float> >&   flow,
	std::vector < std::vector<float> >& CapacityMat,
	std::vector <int>&  visited,
	std::vector <int>&  clusstats0,
	std::vector <int>&  clusstats1,
	std::vector<float>& Wmat);

// this func clears all matrices
void ClearAll(std::vector < std::vector <int> >&   sqlat0,
	std::vector < std::vector <int> >& sqlat1,
	std::vector < std::vector <float> >&   flow,
	std::vector < std::vector<float> >& CapacityMat,
	std::vector <int>&  visited,
	std::vector <int>&  clusstats0,
	std::vector <int>&  clusstats1,
	std::vector<float>& Wmat);
