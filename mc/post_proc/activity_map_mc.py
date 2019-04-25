import lattice as lt
import numpy as np
#import matplotlib.pyplot as plt

reps = 1
L = 40
ntemp=15
phi = 1800
si_s=[1600,1700,1800,1900,2000]

def get_fname(itemp,iter_num,occ,disorder,rep):
	return str(itemp)+'-'+str(iter_num)+'-'+str(occ)+'-'+str(disorder)+'-'+str(rep)+'.csv'


directory = '../data_output/'
	
for si in si_s:
	for itemp in range(ntemp):
		occ_sum = np.zeros((L,L),dtype=int)
		for rep1 in range(1,reps+1):
			#print(rep1)
			source_lattice = lt.Lattice(directory+get_fname(itemp,0,si,phi,rep1))
			for rep2 in range(1,reps+1):
				test_lattice = lt.Lattice(directory+get_fname(itemp,0,si,phi,rep2))
				np.add(occ_sum,np.abs(source_lattice-test_lattice))

		act_map = occ_sum/float(L*L)
		np.savetxt('./activity_maps/'+str(si)+'-'+str(itemp)+'.csv',act_map)

#print(act_map)
#plt.imshow(act_map,vmin=np.min(act_map),vmax=np.max(act_map),interpolation='nearest',cmap='summer')
#plt.colorbar()
#plt.imshow(lt.Lattice(directory+'Si-0_0-1800-0.csv'),alpha=0.2,cmap='PuBu')

#plt.show()
