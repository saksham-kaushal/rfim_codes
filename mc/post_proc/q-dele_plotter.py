import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

L=40
N=L*L

mcresults_fname = '../test1.txt'
data=np.genfromtxt(mcresults_fname,dtype=('f8','f8','f8','f8','f8','f8','f8','f8','i4','U40'),names=('delta','omega','en1','magps1','en2','magps2','q','corr_en','itemp','fname'))
#mcresults_fname = '../../../push-relabel/cloud_reqs/mc/40/mc_VM_40.txt'
#data=np.genfromtxt(mcresults_fname,dtype=('f8','f8','f8','f8','f8','f8','f8','f8','i4'),names=('delta','omega','en1','magps1','en2','magps2','q','corr_en','itemp'))
del_es=data['en2']-data['en1']
qs=data['q']
plt.subplot(2,1,1)
plt.xlim(left=np.min(qs)-0.1,right=1)
plt.plot(qs,del_es,'.')
plt.subplot(2,1,2)
#plt.hist(qs)
plt.xlim(left=np.min(qs)-0.1,right=1)
plt.hist(qs,bins=30)
#sns.distplot(qs,bins=30,hist=True,norm_hist=True)
plt.show()
#plt.savefig('./q-dele_plots/q_dele_plot_'+str(mcresults_fname.split('/')[-1])+'.png')
