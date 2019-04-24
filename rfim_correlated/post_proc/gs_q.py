import numpy as np
import lattice as lt
import os

L=40
N=L*L

directory='../data/'
with open('./gs_q_values/qvals.txt','w') as outfile:
	for fname1 in os.listdir(directory):
		if 'Si' in fname1:
			for fname2 in os.listdir(directory):
				if 'Si' in fname2:
					lt1 = lt.convert_01s(lt.Lattice(directory+fname1))
					lt2 = lt.convert_01s(lt.Lattice(directory+fname2))
					x=np.multiply(lt1,lt2)
					print(len(lt1),len(lt2))
					q=(np.sum(x))/float(N)
					outfile.write(str(fname1.split('.')[0].split('-')[2])+'\t'+str(fname2.split('.')[0].split('-')[2])+'\t'+str(q)+'\t'+str(fname1)+'\t'+str(fname2)+'\n')
