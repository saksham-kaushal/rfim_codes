// g++ -o main LatticeStuff.cpp main_lattice.cpp PushRelabel.cpp -Wall
// ./main

#include "Headers.h"
#include "Parameters.h"
#include "LatticeStuff.h"
#include "PushRelabel.h"

using namespace std;


void core(vector<float>& Bmat, double& maxflow, double& energy, double& mag_ps, vector<int>& visited, string phi_str, string si_str, double delta, int iter, double omega, int rep) 
{
	vector <vector<int>> Exmat(N, vector<int>(N, 0));
	vector < vector <float> > CapacityMat(V, vector<float>(V, 0)), flow(V, vector<float>(V, 0));
	vector<float> Wmat(N, 0);
	
//	cout<<"delta\tGround state energy\tMagnetization per site"<<endl;
	
		
//		create_Bmat_bimodal(Bmat);
		
//		printMatrix(Exmat,N);
//		cout<<endl;
//		cout<<"Disorder "<<iter+1<<endl;
		
			for (auto &i : Exmat)
    			{ fill(i.begin(),i.end(),0); }
			create_Exmat(Exmat);
			fill(Wmat.begin(),Wmat.end(),0);
			for (auto &i : CapacityMat)
    			{ fill(i.begin(),i.end(),0); }
    		for (auto &j : flow)
    			{ fill(j.begin(),j.end(),0); }
			create_CapacityMat(CapacityMat,Exmat);
//			update_Bmat_bimodal(Bmat,delta);
			create_Wmat(Wmat, CapacityMat, Bmat, delta);
			create_Augumented_CapacityMat(Wmat, CapacityMat);
			maxflow = pushRelabel(CapacityMat, flow, 0, V - 1);
//			cout<<"maxflow = "<<maxflow<<endl;
			create_Residual_graph(CapacityMat, flow);
			energy = Hamiltonian(maxflow, CapacityMat);
			dfs(flow,0,visited);
//			savedata(visited,iter,delta,"visited");
			mag_ps = Magnetization(visited);
			
//			savedata(CapacityMat,iter,delta,"Cap");
			savedata(Bmat,iter,delta,phi_str,omega,rep);
			savedata(visited,iter,delta,si_str,omega,rep);
//			savedata(flow,iter,delta,"flow");
		
	
	
	
	
	
	
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
	

	

	Exmat.clear();
	CapacityMat.clear();
	Wmat.clear();
	flow.clear();
}


int main(void)
{
	vector<float>  Bmat1(N,0),Bmat2(N,0);
	vector<int> visited1(V,0),visited2(V,0);
	double delta, maxflow1, energy1, mag_ps1, maxflow2, energy2, energy3, mag_ps2, omega=0.0;
	int iter=0;
	for(iter=0;iter<ITERS;iter++)
	{
		for (delta=delta_beg ; delta<=delta_end; delta+=delta_inc)
		{
			omega = 0.0;
			fill(Bmat1.begin(),Bmat1.end(),0);
			fill(visited1.begin(),visited1.end(),0);
			create_Bmat_gaussian(Bmat1,delta);
			core(Bmat1,maxflow1,energy1,mag_ps1,visited1,"Phi","Si",delta,iter,omega,0);
//			savedata(Bmat1,iter,delta,"Phi",omega);
//			savedata(visited1,iter,delta,"Si",omega);

			for (omega=omega_beg; omega<omega_end; omega+=omega_inc)
			{
				for (int rep=0; rep<9; rep++)
				{
					fill(Bmat2.begin(),Bmat2.end(),0);
					fill(visited2.begin(),visited2.end(),0);
					energy2=0.0;
					energy3=0.0;
					mag_ps2=0.0;
					random_device rand;
					mt19937 generator(rand());
					uniform_real_distribution<> distribution(-1*omega,omega); 
					for (int i=0;i<N;i++)
					{
						Bmat2[i]= Bmat1[i]+distribution(generator);
					}
					core(Bmat2,maxflow2,energy2,mag_ps2,visited2,"Phi_new","Si_new",delta,iter,omega,rep+1);
					energy3 = Hamiltonian_formula(Bmat1,visited2);
					double q = get_q(visited1,visited2);		cout<<fixed_float(delta)<<"\t"<<fixed_float(omega)<<"\t"<<fixed_float(energy1)<<"\t"<<fixed_float(mag_ps1)<<"\t"<<fixed_float(energy2)<<"\t"<<fixed_float(mag_ps2)<<"\t"<<fixed_float(q)<<"\t"<<fixed_float(energy3)<<endl;
				}
			}
		}
	}
	Bmat1.clear();
	Bmat2.clear();
	visited1.clear();
	visited2.clear();
}
