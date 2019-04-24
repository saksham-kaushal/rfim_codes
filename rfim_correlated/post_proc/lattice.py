import numpy as np

def pad_with(vector, pad_width, iaxis, kwargs):
      pad_value = kwargs.get('padder', 0)
      vector[:pad_width[0]] = pad_value
      vector[-pad_width[1]:] = pad_value
      return vector

def convert_01s(occ):
	occ[occ==0] = -1
	return occ
	
def pad(occ,n):
	return np.pad(occ, 1, pad_with,padder=0)

def Lattice(occ_file):
	occ1 = np.genfromtxt(occ_file,delimiter=',')
	occ = np.nan_to_num(occ1)
	occ = np.delete(occ, len(occ), axis=1) 
	return occ
	

