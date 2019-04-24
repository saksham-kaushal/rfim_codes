#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Headers.h"

using namespace std;

void readData(vector<vector<float>> &phi, vector<vector<int>> &si, string filePhi, string fileSi);

double Hamiltonian(vector<vector<float>> &phi, vector<vector<int>> &si);

double Magnetization(vector<vector<int>> &si);
