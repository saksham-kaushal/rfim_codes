#include "Headers.h"
#include "Parameters.h" // >> contains all the req parameters
#include <iostream>
#include <iomanip>


void create_Bmat_bimodal(std::vector <float>& Bmat)
{
	int i = 0, count1 = 0, count2 = 0;
	int var = 0;
//	uint_fast32_t seed = 123456789;
	std::random_device seed_device;
	std::mt19937 generator(seed_device());
	std::uniform_int_distribution<long>  distr(0,10000);
//	std::random_device  rand_dev;
//	std::mt19937        generator(rand_dev());
//	std::cout << "B mtrix ->" << tab;
	for (i = 0; i < N; i++)
	{
		var = distr(generator);
//		std::cout<<var<<" ";
		if (var < 5000)
		{
			Bmat[i] = 1;
//			count1++;
		}
		else
		{
			Bmat[i] = -1;
//			count2++;
		}
//		std::cout << Bmat[i] <<" ";
	}	
//	std::cout << std::endl << "counts : " << count1 << d << count2 << std::endl;


}

//
void update_Bmat_bimodal(std::vector <float>& Bmat, double del)
{
	for (int i = 0; i < N; i++)
	{
		
		if (Bmat[i] > 0)
		{
			Bmat[i] = del;
			
		}
		else
		{
			Bmat[i] = -del;
			
		}
		//std::cout << Bmat[i] << d;
	}
}
// http://www.cplusplus.com/reference/random/normal_distribution/
void create_Bmat_gaussian(std::vector <float>& Bmat,double del)
{
	float var = 0;
	int i = 0,  count1 = 0, count2 = 0;
	std::normal_distribution<double>  distr(0, del);
	std::random_device  rand_dev;
	std::mt19937        generator(rand_dev());
//	std::cout << "B mtrix ->" << tab;
	for (i = 0; i < N; i++)
	{
		//Bmat[i] = distr(generator);// floor(rand() / 1000);
		var = (distr(generator));
		Bmat[i] = var;
//		if (var < 0)
//		{
//			//Bmat[i] = w;
//			count1++;
//		}
//		else
//		{
//			//Bmat[i] = -w;
//			count2++;
//		}
//		std::cout << Bmat[i] << d;
	}	
//	std::cout << std::endl << "counts : " << count1 << d << count2 << std::endl;


}

void create_Exmat(std::vector < std::vector <int> >& Exmat)
{
	int i = 0;
	for (i = 0; i < N; i++)
	{
		if (((i + 1) % VER != 0))
		{
			Exmat[i][i + 1] = 1;
			//cout << Exmat[i][i + 1];
		}
		else
		{
			Exmat[(((i + 1)/VER)-1)*VER][i] = 1;
//			Exmat[i][(((i + 1)/VER)-1)*VER] = 1;
		}
		if (i + VER < N)
		{
			Exmat[i][i + VER] = 1 ;
		}
		else
		{
//			Exmat[i][(i%VER)] = 1;
			Exmat[(i%VER)][i] = 1;
		} 
	}
}


void create_CapacityMat(std::vector < std::vector <float> >& CapacityMat, std::vector < std::vector <int> >& Exmat)
{
	for (int i = 1; i <= N; i++)
	{
		//cout << "\n";
		for (int j = 1; j <= N; j++)
		{
//			if (i < j)
			{
				CapacityMat[i][j] = 4 * J * float(Exmat[i-1][j-1]) ;
			}
/*			else
			{
				CapacityMat[i][j] = 0;
			}
			//cout << CapacityMat[i][j] << " ";
			if (i==VER-1 || j==VER-1)
			{
				CapacityMat[i][j] = 4 * J * float(Exmat[i-1][j-1]) ;
			}
*/		}
	}
}


void create_Wmat(std::vector < float >& Wmat,std::vector < std::vector <float> >& CapacityMat, std::vector < float >& Bmat, double del )
{
	int i = 0, j = 0;
	float cap = 0;
	for (i = 1; i <= N; i++)
	{
		for (j = 1; j <= N; j++)
		{
			cap += CapacityMat[i][j] - CapacityMat[j][i];
			
		}
		//std::cout << cap<<"\t"; 
		Wmat[i - 1] = -2.0 * Bmat[i - 1] - cap / 2.0;
		cap = 0;
		//std::cout <<Wmat[i-1]<< std::endl;
	}
}

// making reisdual graph by overwriting the flow mat,flow[i][j] = CapacityMat[i][j] - flow[i][j];
void create_Augumented_CapacityMat(std::vector < float >& Wmat, std::vector < std::vector <float> >& CapacityMat)
{
	int i = 0;
	for (i = 1; i <= N; i++)
	{
		if (Wmat[i - 1] > 0)
		{
			CapacityMat[0][i] = 0;
			CapacityMat[i][N + 1] = Wmat[i - 1];
		}

		else
		{
			CapacityMat[0][i] = -Wmat[i - 1];
			CapacityMat[i][N + 1] = 0;
		}

	}
}

void create_Residual_graph(std::vector < std::vector <float> >& CapacityMat, std::vector < std::vector <float> >& flow)
{
	int i = 0, j = 0;
	for (i = 0; i < V; i++)
	{
		for (j = 0; j < V; j++)
		{
			flow[i][j] = CapacityMat[i][j] - flow[i][j];
			//cout<< CapacityMat[i][j]-flow[i][j];
		}
		//cout<<endl;
	}
}

double Hamiltonian(double maxflow, std::vector<std::vector<float>>& CapacityMat)
{
	int i,j;
	float corner_cap=0, intermediate_cap=0;
	for (i=1;i<=N;i++)
	{
		corner_cap+=CapacityMat[0][i];
		corner_cap+=CapacityMat[i][N+1];
		for (j=1;j<=N;j++)
		{
			intermediate_cap+=CapacityMat[i][j];
			intermediate_cap+=CapacityMat[j][i];
		}
	}
//	std::cout<<"intermediate cap = "<<intermediate_cap<<std::endl;
//	std::cout<<"corner cap = "<<corner_cap<<std::endl;
//	std::cout<<"maxflow = "<<maxflow<<std::endl;
//	std::cout<<"corner_cap "<<corner_cap<<" intermediate_cap "<<intermediate_cap<<std::endl;
	return (maxflow-(corner_cap/2.0)-(intermediate_cap/8.0));
}

double Magnetization(std::vector <int>& visited)
{
	long mag = 0;
	double mag_ps = 0.0;
	for(int i=1;i<=N;i++)
	{
		if (visited[i]==0)
		{
			mag-=1;
		}
		else if (visited[i]==1)
		{
			mag+=1;
		}
	}
	mag_ps = double(mag)/double(N);
	return mag_ps;
}

void printMatrix(std::vector <std::vector <int> >  & M, int len) {
	int i, j;
	for (i = 0; i < len; i++) {
		for (j = 0; j < len; j++)
			std::cout << " " << std::setw(3) << M[i][j]; 
		std::cout << "\n";
	}
}

void printMatrix(std::vector <std::vector <float> >  & M, int len) {
	int i, j;
	for (i = 0; i < len; i++) {
		for (j = 0; j < len; j++)
			std::cout << " " << std::setw(3) << M[i][j];
		std::cout << "\n";
	}
}

void printMatrix(std::vector<float> & M, int len) {
	int i;
	for (i = 0; i < len; i++) 
	{
		std::cout << " " << std::setw(3) << M[i];
	}
	std::cout << "\n";
}

void printMatrix(std::vector<int> & M, int len) {
	int i;
	for (i = 0; i < len; i++) 
	{
		std::cout << " " << std::setw(3) << M[i];
	}
	std::cout << "\n";
}


void savedata(std::vector< std::vector< float > > & Mat, int l, double del, std::string s)
{
	int delta1 = int(del * 1000);
	
//	std::random_device rand_seed;
//	std::mt19937 generator(rand_seed());
//	std::uniform_int_distribution<long>  distribution(0,100);
//	long rand_num = distribution(generator);
	
	std::ofstream fobj("./data/"+s + "-" + std::to_string(l) + "-" + std::to_string(delta1) + ".txt");
//	std::ofstream fobj(s + "-" + std::to_string(l) + "-" + std::to_string(delta1) + "_" + std::to_string(rand_num) + ".txt");
	//fobj << "iteration : " << l << "\ndelta : " << del << std::endl;
	for (int i = 0; i < int(Mat.size()); i++)
	{
		for (int j = 0; j < int(Mat[i].size()); j++)
		{
//			if (Mat[i][j] != 0)
			{
//				fobj <<i<<" "<<j<<" "<< Mat[i][j]<<" \n";
				fobj << Mat[i][j] <<" ";
			}
			
		}
		fobj << std::endl;
	}
	fobj.clear();
}

void savedata(std::vector< int > & Mat,int l,double del,std::string s)
{
	int delta1 = int(del * 1000);
	
//	std::random_device rand_seed;
//	std::mt19937 generator(rand_seed());
//	std::uniform_int_distribution<long>  distribution(0,100);
//	long rand_num = distribution(generator);
	
	std::ofstream fobj("./data/"+s+"-" + std::to_string(l) + "-" + std::to_string(delta1) + ".csv");
//	std::ofstream fobj(s+"-" + std::to_string(l) + "-" + std::to_string(delta1) + "_" + std::to_string(rand_num) + ".csv");
	//fobj << "iteration : " << l << "\ndelta : " << del << std::endl;

//	Mat.erase(Mat.begin() + 1); 		//Probable reason of error
//	Mat.erase(Mat.begin());

	for (int i = 1; i <= N; i++)
	{
		if (i%VER == 1 && i!=1) fobj << std::endl;
		fobj << Mat[i] << ",";
	}
	fobj.clear();
}

void savedata(std::vector< float > & Mat, int l, double del, std::string s)
{
	int delta1 = int(del * 1000);
	
//	std::random_device rand_seed;
//	std::mt19937 generator(rand_seed());
//	std::uniform_int_distribution<long>  distribution(0,100);
//	long rand_num = distribution(generator);
	
	std::ofstream fobj("./data/"+s + "-" + std::to_string(l) + "-" + std::to_string(delta1) + ".csv");
//	std::ofstream fobj(s + "-" + std::to_string(l) + "-" + std::to_string(delta1) + "_" + std::to_string(rand_num) + ".csv");
	//fobj << "iteration : " << l << "\ndelta : " << del << std::endl;
	
	
		
		for (int i = 0; i < N; i++)
		{

			if (i%VER == 0 && i)	fobj << std::endl;
			fobj << Mat[i] << ",";

		}
	

	fobj.clear();
}



//-----------------------------------------------------------------

void dfs(std::vector <std::vector <float> >& F, int s, std::vector <int>& visited)
{
	visited[s] = 1;
	for (int i = 0; i < V; i++)
	{
		if (F[s][i]>0 && !visited[i])
		{
			//std::cout << i << d;
			dfs(F, i, visited);
		}
	}

}

void getSi(std::vector <int>& visited, int len)
{
	int i=0;
	for (i=1;i<len-1;i++)
	{
		std::cout << std::setw(3) << 2*visited[i]-1;
		if ((i)%VER == 0) {std::cout<<std::endl;}
	}
}
