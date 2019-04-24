#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include "Parameters.h"
#include "HamiltonianTest.h"

using namespace std;

int main(void)
{
	double delta;
	int i;
	for(i=0;i<ITERS;i++)
	{
		for(delta=delta_beg;delta<=delta_end;delta+=delta_inc)
		{
			long num = static_cast<long>(1000*delta);
			string phiFile="./Phi-"+to_string(i)+"-"+to_string(num)+".csv";
			string siFile="./Si-"+to_string(i)+"-"+to_string(num)+".csv";
			double tot_en=0.0;
			vector<vector<float>> phi(VER,vector<float>(VER,0));
			vector<vector<int>> si(VER,vector<int>(VER,0));
			readData(phi,si,phiFile,siFile);
			tot_en=Hamiltonian(phi,si);
			double mag_per_site=Magnetization(si);
			cout<<phiFile<<siFile<<endl;
			cout<<"total energy and magnetization per site for delta = "<<delta<<" are "<<tot_en<<" and "<<mag_per_site<<endl<<endl;
		}
	
	}
	return 0;
}

void readData(vector<vector<float>> &phi, vector<vector<int>> &si, string filePhi, string fileSi)
{
	ifstream infilePhi, infileSi;
	infilePhi.open(filePhi);
	infileSi.open(fileSi);
	if(!infilePhi || !infileSi)
	{
		cerr<<"cannot open file "<<filePhi<<", "<<fileSi<<endl;
//		exit(1);
	}
	float phis;
	int sis;
	char coma;
	int i,j;
	for (i=0;i<VER;i++)
	{
		for(j=0;j<VER;j++)
		{
			infilePhi>>phis>>coma;
			infileSi>>sis>>coma;
			phi[i][j]=phis;
			if(sis<1)
				{ si[i][j]=-1; }
			else
				{ si[i][j]=1; }
		}
	}
	infilePhi.close();
	infileSi.close();
}


double Hamiltonian(vector<vector<float>> &phi, vector<vector<int>> &si)
{
	double RF_component=0.0;
	long IM_component=0;
	unsigned int i,j;
	for(i=0;i<VER;i++)
	{
		for(j=0;j<VER;j++)
		{
			RF_component+=(phi[i][j]*(si[i][j]));
			IM_component+=(si[i][j]*(si[(i-1)%VER][j]+si[i][(j-1)%VER]+si[(i+1)%VER][j]+si[i][(j+1)%VER]));
		}
	}
	//si[i-1][j]+si[i][j-1]+si[i+1][j]+si[i][j+1]
	cout<<"RF_component "<<RF_component<<endl<<"IM_component "<<IM_component<<endl;
	return (-(J*(double(IM_component))/2.0)-RF_component);
}


double Magnetization(vector<vector<int>> &si)
{
	int i,j;
	double mag,mag_per_site;
	for(i=0;i<VER;i++)
	{
		for(j=0;j<VER;j++)
		{
			mag+=si[i][j];
		}
	}
//	cout<<"Net magnetization = "<<mag<<endl;
	mag_per_site = mag/double(VER*VER);
	return mag_per_site;
}
