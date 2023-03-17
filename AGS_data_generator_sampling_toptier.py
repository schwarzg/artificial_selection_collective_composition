import numpy as np
from growth_alter import *
from selection import select_toptier
from selection import community_function as cf
import os,sys

#import matplotlib.pyplot as plt

#Model parameter
mu=1e-4
r=0.5
s=3e-2
N0=1000
mbar=50
ncomm=100
rhat=0.1
nensemble=3
ncycle=1000
tcycle=np.log(ncomm+1)/r
nsel=5
nrep=int(ncomm/nsel)

if len(sys.argv)!=1:
	mbar=int(sys.argv[1])	
	rhat=float(sys.argv[2])

descriptor="AGS_PD_samp_N0%s_mbar%s_r%s_s%s_mu%s_ncomm%d_nsel%d_rhat%s_ncycle%d/"%(N0,mbar,r,s,mu,ncomm,nsel,rhat,ncycle) 
folder="data/raw/"+descriptor

nproc=1
rank=0

'''
#MPI start when mpi4py is imported
from mpi4py import MPI

comm = MPI.COMM_WORLD
nproc = comm.Get_size()
rank = comm.Get_rank()
'''

if rank==0:
	m0sel_ini=np.sort(np.clip(np.random.binomial(N0,mbar/N0,size=ncomm).astype(int),0,N0))
else:
	m0sel_ini=np.empty(ncomm,dtype=int)

#remove comment if use mpi
#m0sel_ini=comm.bcast(m0sel_ini,root=0)

w0sel_ini=N0-m0sel_ini
	
print("Initial setting:",w0sel_ini,m0sel_ini)

if rank==0:
	if not os.path.exists(folder):
		os.mkdir(folder)
#comm.barrier()

for e in range(rank,nensemble,nproc):
	#prepare initial state
	w0sel=np.copy(w0sel_ini)
	m0sel=np.copy(m0sel_ini)
	
	#prepare data space
	t_selected=np.array([])
	w_selected=np.array([])
	m_selected=np.array([])
	ind_selected=np.array([])
	
	#Artificial selection cycle : Growth, Selection, and Reproduction
	#run community lifecycle
	for i in range(ncycle):
		w_selected=np.append(w_selected,w0sel)
		m_selected=np.append(m_selected,m0sel)
	
		#Space for selection	
		lastw=np.zeros(ncomm)
		lastm=np.zeros(ncomm)

		#For each community
		for j in range(ncomm):
			#store grown number for selection
			lastw[j]=growth_sampling_w(w0sel[j],0,tcycle,r,mu,s)
			lastm[j]=growth_sampling_m(m0sel[j],0,w0sel[j],0,0,tcycle,r,mu,s)
	
		#Selection phase
		w_selected=np.append(w_selected,lastw)
		m_selected=np.append(m_selected,lastm)
		ind_sel=select_toptier(nsel,lastw,lastm,r,s,rhat)
		w_sel=lastw[ind_sel]
		m_sel=lastm[ind_sel]
		ind_selected=np.append(ind_selected,ind_sel)
		t_selected=np.append(t_selected,(i+1)*tcycle)	
	
		#reproduction phase
		for i,ind in enumerate(ind_sel):	
			m0sel[nrep*i:nrep*(i+1)]=np.random.binomial(N0,cf(w_sel[i],m_sel[i]),size=nrep)
		m0sel=np.sort(m0sel)
		w0sel=N0-m0sel
	
	#Save data
	#np.savetxt for data at cycle time only
	output=np.hstack((np.array([t_selected]).T,ind_selected.reshape(ncycle,nsel),w_selected.reshape(ncycle,2*ncomm),m_selected.reshape(ncycle,2*ncomm)))
	np.savetxt(folder+"%d.cycle"%(e),output)

#Finish MPI
#MPI.COMM_WORLD.Barrier()
#MPI.Finalize()


