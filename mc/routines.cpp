#include "headers.h"
#include "parameters.h"

void readSi(std::vector<std::vector<int>>& occ, double si, int iter)
{
	std::string fname = "./data_input/Si-0_"+std::to_string(iter)+"-"+std::to_string(int(si*10+0.5)*100)+".csv";
	std::ifstream infile;
	infile.open(fname);
	if(!infile)
	{
		std::cerr<<"cannot open file "<<fname<<std::endl;
//		exit(1);
	}
	int config;
	char coma;
	for(int i=0;i<VER;i++)
	{
		for(int j=0;j<VER;j++)
		{
			infile>>config>>coma;
			if(config<1)
				{ occ[i][j]=-1; }
			else
				{ occ[i][j]=1; }
		}
	}
	infile.close();
}

void readPhi(std::vector<std::vector<float>>& phi, float si, int iter)
{
	std::string fname = "./data_input/Phi-0-"+std::to_string(iter)+"-"+std::to_string(int(si*10+0.5)*100)+".csv";
	std::ifstream infile;
	infile.open(fname);
	if(!infile)
	{
		std::cerr<<"cannot open file "<<fname<<std::endl;
//		exit(1);
	}
	float phis;
	char coma;
	for(int i=0;i<VER;i++)
	{
		for(int j=0;j<VER;j++)
		{
			infile>>phis>>coma;
			phi[i][j]=phis;
		}
	}
	infile.close();
}



double get_q(std::vector<std::vector<int>>& occ_init, std::vector<std::vector<int>>& occ)
{
	long q = 0;
	long x = 0;
	for(int i=0;i<VER;i++)
	{
		for(int j=0; j<VER; j++)
		{
			x=(occ_init[i][j])*(occ[i][j]);
			q+=x;
		}
	}
	return (double(q)/(VER*VER));
}


double Hamiltonian(std::vector<std::vector<int>>& si, std::vector<std::vector<float>>& phi)
{
	double RF_component=0.0;
	long IM_component=0;
	int i,j;
	for(i=0;i<VER;i++)
	{
		for(j=0;j<VER;j++)
		{
			int pre_i, pre_j, post_i, post_j;
			if (i==0) {pre_i = VER-1;}
			else {pre_i=i-1;}
		
			if (j==0) {pre_j = VER-1;}
			else {pre_j=j-1;}
			
			if (i==VER-1) { post_i = 0; }
			else {post_i=i+1;}
	
			if (j==VER-1) { post_j = 0; }
			else {post_j=j+1;}
			
			RF_component+=(phi[i][j]*(si[i][j]));
			IM_component+=(si[i][j]*(si[pre_i][j]+si[i][pre_j]+si[post_i][j]+si[i][post_j]));
		}
	}
	//si[i-1][j]+si[i][j-1]+si[i+1][j]+si[i][j+1]
//	cout<<"RF_component "<<RF_component<<endl<<"IM_component "<<IM_component<<endl;
	return (-(J*(double(IM_component))/2.0)-RF_component);
}


double Magnetization(std::vector<std::vector<int>> &si)
{
	int i,j;
	double mag=0,mag_per_site=0;
	for(i=0;i<VER;i++)
	{
		for(j=0;j<VER;j++)
		{
			mag+=si[i][j];
		}
	}
//	cout<<"Net magnetization = "<<mag<<endl;
//	mag_per_site = mag/double(VER*VER);
	return mag;
}


void readBeta(std::vector<float>& beta)
{
	std::string fname = "./data_input/beta.txt";
	std::ifstream infile;
	infile.open(fname);
	if(!infile)
	{
		std::cerr<<"cannot open file "<<fname<<std::endl;
//		exit(1);
	}
//	float p;
	for(int i=0;i<NTEMP;i++)
	{
		float betaa;
		infile>>betaa;
		beta[i]=betaa;
	}
	infile.close();
}


double del_h(std::vector<std::vector<int>>& si, std::vector<std::vector<float>>& phi, int ic, int jc)
{
	int pre_ic, pre_jc, post_ic, post_jc;
	if (ic==0) {pre_ic = VER-1;}
	else {pre_ic=ic-1;}
	
	if (jc==0) {pre_jc = VER-1;}
	else {pre_jc=jc-1;}
	
	if (ic==VER-1) { post_ic = 0; }
	else {post_ic=ic+1;}
	
	if (jc==VER-1) { post_jc = 0; }
	else {post_jc=jc+1;}
//	std::cout<<ic<<" "<<jc<<" "<<pre_ic<<" "<<pre_jc<<" "<<post_ic<<" "<<post_jc<<std::endl;
	return 2*(si[ic][jc]*phi[ic][jc])+2*si[ic][jc]*(si[pre_ic][jc]+si[ic][pre_jc]+si[post_ic][jc]+si[ic][post_jc]);
}

double del_m(std::vector<std::vector<int>>& si, int ic, int jc)
{
	return 2*si[ic][jc];
}

double update_h(double toten,double delta_h)
{
	toten+=delta_h;
	return toten;
}

double update_mag(double mag,double delta_m)
{
	mag-=delta_m;
	return mag;
}

void savedata(std::vector<std::vector<int>> & Mat,int itemp, float si, float phi,int iter, int rep)
{
	long delta1 = long(si * 1000);
	long delta2 = long(phi * 1000);
//	std::random_device rand_seed;
//	std::mt19937 generator(rand_seed());
//	std::uniform_int_distribution<long>  distribution(0,100);
//	long rand_num = distribution(generator);
	
	std::string fname = "./data_output/" + std::to_string(itemp) + "-" + std::to_string(iter) + "-" + std::to_string(delta1) + "-" + std::to_string(delta2) + "-" + std::to_string(rep) + ".csv";
	std::cout<<fname<<"\n";
	std::ofstream fobj(fname);
//	std::ofstream fobj(s+"-" + std::to_string(l) + "-" + std::to_string(delta1) + "_" + std::to_string(rand_num) + ".csv");
	//fobj << "iteration : " << l << "\ndelta : " << del << std::endl;

	for(int i=0;i<VER;i++)
	{
		for(int j=0;j<VER;j++)
		{
			if(Mat[i][j]<0)
			{
				fobj<<"0,";
			}
			else
			{
				fobj<<"1,";
			}
			
		}
		fobj<<std::endl;
	}
	fobj.clear();
}

