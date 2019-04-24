#include "headers.h"
#include "routines.h"
#include "parameters.h"

int main(void)
{
	std::vector<float> beta(NTEMP);
	readBeta(beta);
	
	for (float disorder=delta_beg; disorder<delta_end; disorder+=delta_inc)
	{
		for (int iter=0;iter<ITERS;iter++)
		{
			std::vector<std::vector<int>> occ_init(VER,std::vector<int>(VER,0));
			readSi(occ_init,disorder,iter);
			std::vector<std::vector<float>> phi(VER,std::vector<float>(VER,0));
			readPhi(phi,disorder,iter);
			
			double toten_init = Hamiltonian(occ_init,phi);
			double maginit = Magnetization(occ_init);
			double magps_init = maginit/float(NS);
			//std::cout<<maginit<<magps_init<<std::endl;
			
			std::random_device seed_device;
			std::mt19937 generator_real(seed_device());
			std::uniform_real_distribution<float> distribution_real(0,1);
			std::mt19937 generator_int(seed_device());
			std::uniform_int_distribution<int> distribution_int(0,VER-1);
			
			for (float si=si_beg; si<si_end; si+=si_inc)
			{
				//for (int iter2=0;iter2<ITERS;iter2++)
				{
					std::vector<std::vector<int>> occ(VER,std::vector<int>(VER,0));
					readSi(occ,si,iter);
				
					toten = Hamiltonian(occ,phi);
					mag = Magnetization(occ);
					
					for (int rep=0; rep<1000; rep++)
					{
					double delta_h=0.0, delta_m=0.0, toten=0.0, mag=0.0, q=0.0;
					
					for(int itemp=0;itemp<NTEMP;itemp++)
					{
						for(int imeas=0;imeas<int((NMEAS*100)/beta[itemp]);imeas++)
						{
							for(int nmcs=0;nmcs<NS;nmcs++)
							{
								int ic = distribution_int(generator_int);
								int jc = distribution_int(generator_int);
								//std::cout<<ic<<jc<<std::endl;
								unsigned int ix = ic;
								unsigned int jx = jc;
								delta_h = del_h(occ,phi,ix,jx);
								delta_m = del_m(occ,ic,jc);
								if(delta_h < 0.0)
								{
									occ[ic][jc] = -occ[ic][jc];
									toten = update_h(toten,delta_h);
									mag = update_mag(mag,delta_m);
								}
								else
								{
									float a = exp(-delta_h*beta[itemp]);
									float c = distribution_real(generator_real);
									if(c<a)
									{
										occ[ic][jc] = -occ[ic][jc];
										toten = update_h(toten,delta_h);
										mag = update_mag(mag,delta_m);
									}
								}	
							}	//nmcs loop
						}	//imeas loop
						float magps = mag/float(NS);
						q = get_q(occ_init,occ);
						std::cout<<fixed_float(disorder)<<"\t"<<fixed_float(beta[itemp])<<"\t"<<fixed_float(toten_init)<<"\t"<<fixed_float(magps_init)<<"\t"<<fixed_float(toten)<<"\t"<<fixed_float(magps)<<"\t"<<fixed_float(q)<<"\t"<<fixed_float(si)<<"\t"<<fixed_float(itemp)<<std::endl;
						savedata(occ,itemp,si,disorder,iter,rep+1);
					}	//rep loop
					}	//itemp loop
			occ.clear();
				}	//iter2 loop
			}	//si loop
		phi.clear();
		occ_init.clear();
		}	//iter loop
	beta.clear();
	}	//delta loop
}
