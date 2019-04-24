#pragma once

#include "headers.h"

void readSi(std::vector<std::vector<int>>& occ, double si, int iter);
void readPhi(std::vector<std::vector<float>>& phi, float disorder, int iter);
double get_q(std::vector<std::vector<int>>& occ_init, std::vector<std::vector<int>>& occ);
double Hamiltonian(std::vector<std::vector<int>>& occ, std::vector<std::vector<float>>& phi);
double Magnetization(std::vector<std::vector<int>> &si);
void readBeta(std::vector<float>& beta);
double del_h(std::vector<std::vector<int>>& si, std::vector<std::vector<float>>& phi, int ic, int jc);
double del_m(std::vector<std::vector<int>>& si, int ic, int jc);
double update_h(double toten,double delta_h);
double update_mag(double mag,double delta_m);
void savedata(std::vector<std::vector<int>> & Mat,int itemp, float si, float phi, int iter, int rep);
