import numpy as np
from model import model
from tau_leaping import tau_leaping
from selection import select
from selection import community_function as cf
from selection import score_function as sf
from reproduction import hypergeometric_reproduce
import h5py as h5
import os,sys

import matplotlib.pyplot as plt
#initial setting
print("Arguement : mbar")
mu=1e-4
r=0.5
s=5e-3
N0=1000
mbar=100
ncomm=10
rhat=0.05

if len(sys.argv)!=1:
	mbar=int(sys.argv[1])	

proFunc=[
lambda x: r*x[0],
lambda x: mu*x[0],
lambda x: (r+s)*x[1]
]
changeVec=np.array([[1,0],[-1,1],[0,1]]).astype(int)
rxnOrder=np.array([[1,0],[1,0],[0,1]]).astype(int)



nensemble=100
ncycle=30
tcycle=np.log(ncomm+1)/r

#MPI start when mpi4py is imported
from mpi4py import MPI

comm = MPI.COMM_WORLD
nproc = comm.Get_size()
rank = comm.Get_rank()

if rank==0:
	m0sel_ini=np.clip(np.random.normal(loc=mbar,scale=np.sqrt(mbar),size=ncomm).astype(int),0,1000)
else:
	m0sel_ini=np.empty(ncomm,dtype=int)

m0sel_ini=comm.bcast(m0sel_ini,root=0)

w0sel_ini=N0-m0sel_ini
	
print("Initial setting:",w0sel_ini,m0sel_ini)
growthmodel=model(proFunc=proFunc,changeVec=changeVec,rxnOrder=rxnOrder)

solver=tau_leaping(model=growthmodel,eps=0.01)

w0sel=np.copy(w0sel_ini)
m0sel=np.copy(m0sel_ini)


folder="data/raw/"

if rank==0:
	if not os.path.exists(folder):
		os.mkdir(folder)
comm.barrier()

for e in range(rank,nensemble,nproc):
	
	#prepare data space
	tsel=np.array([])
	wsel=np.array([])
	msel=np.array([])
	
	t_selected=np.array([])
	w_selected=np.array([])
	m_selected=np.array([])
	
	#GR model : Growth and Reproduction
	#run community lifecycle
	for i in range(ncycle):
	
		#Space for selection	
		lastw=np.array([])
		lastm=np.array([])

		#For each community
		for j in range(ncomm):

			#Growth phase
			T,X=solver.run(np.array([w0sel[j],m0sel[j]]),0,tcycle)
			tsel=np.append(tsel,T[:-1]+i*tcycle)	
			wsel=np.append(wsel,X[0,:-1])
			msel=np.append(msel,X[1,:-1])

			#store for selection
			lastw=np.append(lastw,X[0,-1])
			lastm=np.append(lastm,X[1,-1])	
			#reproduction phase	
			w0sel[j],m0sel[j]=hypergeometric_reproduce(X[0,-1],X[1,-1],N0,Ngroup=1)	

		#Selection phase
		w_selected=np.append(w_selected,lastw)
		m_selected=np.append(m_selected,lastm)
		t_selected=np.append(t_selected,(i+1)*tcycle)	
	
	
	#Save data
	#np.savetxt for data at cycle time only
	output=np.vstack((np.array([t_selected]),w_selected.reshape(ncomm,ncycle),m_selected.reshape(ncomm,ncycle))).T
	np.savetxt(folder+"%s_%s_%s_%s_%s_%s_%s_NS_%d.cycle"%(N0,mbar,r,s,mu,ncomm,rhat,e),output)

	'''
	if e==0:
		#Prepare datafile for ensemble
		f=h5.File(folder+"%s_%s_%s_%s_%s_%s_%s_NS_%d.hdf5"%(N0,mbar,r,s,mu,ncomm,rhat,e),'w')
	
		#data structure /comm/T,w,m and /s/T,w,m
		for n in range(ncomm):
			com=f.create_dataset('%d'%n,data=np.array([tsel,wsel,msel]).T)
	
		f.close()
	'''
#Finish MPI
MPI.COMM_WORLD.Barrier()
MPI.Finalize()


