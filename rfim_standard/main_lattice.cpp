// g++ -o main LatticeStuff.cpp main_lattice.cpp PushRelabel.cpp -Wall
// ./main

#include "Headers.h"
#include "Parameters.h"
#include "LatticeStuff.h"
#include "PushRelabel.h"

using namespace std;

int main(void) 
{
	vector <vector<int>> Exmat(N, vector<int>(N, 0));
	vector < vector <float> > CapacityMat(V, vector<float>(V, 0)), flow(V, vector<float>(V, 0));
	vector<float>  Bmat(N, 0), Wmat(N, 0);
	vector <int>  visited(V, 0);
	double delta, maxflow, energy, mag_ps;
	int iter=0;
	
//	cout<<"delta\tGround state energy\tMagnetization per site"<<endl;
	for(iter=0;iter<ITERS;iter++)
	{
		
//		create_Bmat_bimodal(Bmat);
		
//		printMatrix(Exmat,N);
//		cout<<endl;
//		cout<<"Disorder "<<iter+1<<endl;
		for (delta=delta_beg ; delta<=delta_end; delta+=delta_inc)
		{
			fill(Bmat.begin(),Bmat.end(),0);
			for (auto &i : Exmat)
    			{ fill(i.begin(),i.end(),0); }
			create_Exmat(Exmat);
			fill(visited.begin(),visited.end(),0);
			fill(Wmat.begin(),Wmat.end(),0);
			for (auto &i : CapacityMat)
    			{ fill(i.begin(),i.end(),0); }
    		for (auto &j : flow)
    			{ fill(j.begin(),j.end(),0); }
			create_CapacityMat(CapacityMat,Exmat);
//			update_Bmat_bimodal(Bmat,delta);
			create_Bmat_gaussian(Bmat,delta);
			create_Wmat(Wmat, CapacityMat, Bmat, delta);
			create_Augumented_CapacityMat(Wmat, CapacityMat);
			maxflow = pushRelabel(CapacityMat, flow, 0, V - 1);
//			cout<<"maxflow = "<<maxflow<<endl;
			create_Residual_graph(CapacityMat, flow);
			energy = Hamiltonian(maxflow, CapacityMat);
			dfs(flow,0,visited);
//			savedata(visited,iter,delta,"visited");
			mag_ps = Magnetization(visited);
			cout<<fixed_float(delta)<<"\t"<<fixed_float(energy)<<"\t\t"<<fixed_float(mag_ps)<<endl;
			
//			savedata(CapacityMat,iter,delta,"Cap");
			savedata(Bmat,iter,delta,"Phi");
			savedata(visited,iter,delta,"Si");
//			savedata(flow,iter,delta,"flow");
		}
	}
	
	
	
	
	
//	printMatrix(Exmat,N);
//	cout<<endl;
//	printMatrix(Bmat,N);
//	cout<<endl;
//	printMatrix(Wmat,N);
//	cout<<endl;
//	printMatrix(CapacityMat,V);
//	printMatrix(CapacityMat,V);
//	cout<<endl;
//	printMatrix(flow,V);
//	cout<<endl;
	
//	printMatrix(flow,V);
//	cout<<endl;
	
//	printMatrix(visited,V);
//	getSi(visited,V);
	

	
	Bmat.clear();
	Exmat.clear();
	CapacityMat.clear();
	Wmat.clear();
	flow.clear();
	visited.clear();
	return 0;
}
