import os
import subprocess as sp

outfiles={}
processes=[]
for i in range(14):
	directory = './'+str(i+1)
	filename = '/randex_80_1800_500-'+str(i+1)+'.txt'
	os.mkdir(directory)
	os.mkdir(directory+'/data')
	sp.call(['cp','main',directory])
	sp.call(['touch',directory+filename])
	outfiles[filename]=open(directory+filename,'w')
	processes.append(sp.Popen(['./main'],stdout=outfiles[filename],cwd=directory))
	
exit_status = [p.wait() for p in processes]
