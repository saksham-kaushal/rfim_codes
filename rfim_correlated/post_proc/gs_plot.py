import lattice as lt
import matplotlib.pyplot as plt
import os

directory = '../data/'
for fname in os.listdir(directory):
	if 'Si' in fname:
		plt.imshow(lt.Lattice(directory+fname))
		plt.title(fname)
		plt.savefig('gs_plots/'+fname+'.png')
		plt.close('all')
