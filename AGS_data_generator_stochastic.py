import numpy as np
from model import model
from tau_leaping import tau_leaping
from selection import select
from selection import community_function as cf
import os,sys


#Model parameter
mu=1e-4
r=0.5
s=5e-2
N0=1000
mbar=100
ncomm=10
rhat=0.05
nensemble=1
ncycle=10
tcycle=np.log(ncomm+1)/r
folder="data/raw/"



if len(sys.argv)!=1:
	mbar=int(sys.argv[1])	
	rhat=float(sys.argv[2])

descriptor="AGS_PD_sto_N0%s_mbar%s_r%s_s%s_mu%s_ncomm%d_rhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,rhat,ncycle) 
folder="data/raw/"+descriptor

proFunc=[
lambda x: r*x[0],
lambda x: mu*x[0],
lambda x: (r+s)*x[1]
]
changeVec=np.array([[1,0],[-1,1],[0,1]]).astype(int)
rxnOrder=np.array([[1,0],[1,0],[0,1]]).astype(int)




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
growthmodel=model(proFunc=proFunc,changeVec=changeVec,rxnOrder=rxnOrder)

solver=tau_leaping(model=growthmodel,eps=0.01)

w0sel=np.copy(w0sel_ini)
m0sel=np.copy(m0sel_ini)



if rank==0:
	if not os.path.exists(folder):
		os.mkdir(folder)
#comm.barrier()

for e in range(rank,nensemble,nproc):
	
	#prepare data space
	tsel=np.array([])
	wsel=np.array([])
	msel=np.array([])
	
	t_selected=np.array([])
	w_selected=np.array([])
	m_selected=np.array([])
	ind_selected=np.array([])
	
	#GSR model : Growth, Selection, and Reproduction
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
	
		#Selection phase
		w_selected=np.append(w_selected,lastw)
		m_selected=np.append(m_selected,lastm)
		ind_sel=select(lastw,lastm,r=r,s=s,rhat=rhat)
		w_sel=lastw[ind_sel]
		m_sel=lastm[ind_sel]
		ind_selected=np.append(ind_selected,ind_sel)
		t_selected=np.append(t_selected,(i+1)*tcycle)	
	
		#reproduction phase	
		m0sel=np.random.binomial(N0,cf(w_sel,m_sel),size=ncomm)
		m0sel=np.sort(m0sel)
		w0sel=N0-m0sel
	
	if e==0:
		plt.show()	
	#Save data
	#np.savetxt for data at cycle time only
	output=np.hstack((np.array([t_selected,ind_selected]).T,w_selected.reshape(ncycle,2*ncomm),m_selected.reshape(ncycle,2*ncomm)))
	np.savetxt(folder+"%d.cycle"%(e),output)
	

	

#Finish MPI
#MPI.COMM_WORLD.Barrier()
#MPI.Finalize()


