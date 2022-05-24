import numpy as np
from model import model
from tau_leaping import tau_leaping
from selection import select
from selection import community_function as cf
from selection import score_function as sf
from reproduction import hypergeometric_reproduce
#import h5py as h5
import os,sys

import matplotlib.pyplot as plt

#Model parameter
print("Arguement : mbar")
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

cset=[

'#a50026',
'#d73027',
'#f46d43',
'#fdae61',
'#fee090',
'#e0f3f8',
'#abd9e9',
'#74add1',
'#4575b4',
'#313695'
]


if len(sys.argv)!=1:
	mbar=int(sys.argv[1])	

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
if rank==0:
	m0sel_ini=np.clip(np.random.normal(loc=mbar,scale=np.sqrt(mbar),size=ncomm).astype(int),0,1000)
else:
	m0sel_ini=np.empty(ncomm,dtype=int)

m0sel_ini=comm.bcast(m0sel_ini,root=0)
'''

m0sel_ini=np.sort(np.clip(np.random.normal(loc=mbar,scale=np.sqrt(mbar),size=ncomm).astype(int),0,1000))
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
			print('cycle',j,'from',T[0],'to',T[-1])
			wsel=np.append(wsel,X[0,:-1])
			msel=np.append(msel,X[1,:-1])

			#store for selection
			lastw=np.append(lastw,X[0,-1])
			lastm=np.append(lastm,X[1,-1])	
	
			#Trajectory data temporary test
			if e==0:
				import matplotlib.pyplot as plt
				cfs=cf(X[0],X[1])
				plt.plot(T+i*tcycle,cfs,c=cset[j])

		#Selection phase
		w_selected=np.append(w_selected,lastw)
		m_selected=np.append(m_selected,lastm)
		ind_sel=select(lastw,lastm,r=r,s=s,rhat=rhat)
		w_sel=lastw[ind_sel]
		m_sel=lastm[ind_sel]
		ind_selected=np.append(ind_selected,ind_sel)
		t_selected=np.append(t_selected,(i+1)*tcycle)	
	
		#reproduction phase	
		w0sel,m0sel=hypergeometric_reproduce(w_sel,m_sel,N0,Ngroup=10)	
		m0sel=np.sort(m0sel)
		w0sel=N0-m0sel
	
	if e==0:
		plt.show()	
	#Save data
	#np.savetxt for data at cycle time only
	output=np.hstack((np.array([t_selected,ind_selected]).T,w_selected.reshape(ncycle,ncomm),m_selected.reshape(ncomm,ncycle)))
	np.savetxt(folder+"%s_%s_%s_%s_%s_%s_%s_AGS_%d.cycle"%(N0,mbar,r,s,mu,ncomm,rhat,e),output)

	#Trajectory data
	#np.savetxt(folder+"%s_%s_%s_%s_%s_%s_%s_AGS_%d.trajectory"%(N0,mbar,r,s,mu,ncomm,rhat,e),np.array([tsel,wsel,msel]))
	

	
	'''
	#Prepare datafile for ensemble
	f=h5.File(folder+"%s_%s_%s_%s_%s_%s_%s_AGS_%d.hdf5"%(N0,mbar,r,s,mu,ncomm,rhat,e),'w')

	#data structure /comm/T,w,m and /s/T,w,m
	for n in range(ncomm):
		com=f.create_dataset('%d'%n,data=np.array([tsel,wsel,msel]).T)
	#sel=f.create_dataset('s',data=np.array([t_selected,w_selected,m_selected]).T)

	f.close()

	import matplotlib.pyplot as plt
	cx2.hlines(rhat,0,ncycle*tcycle,linestyles=':')
	averaged=np.w_selected
	cx2.plot(t_select,averaged,marker='x',ms=5,c='red',ls='--',label='w/ sel., avr')
	cx2.fill_between(tselect,averaged+deviated,averaged-deviated,color='red',alpha=0.2)
	cx2.fill_between(tselect,averagedGR+deviatedGR,averagedGR-deviatedGR,color='blue',alpha=0.2)
	cx2.set_yscale('log')
	cx2.set_ylabel(r'$C$')
	cx2.set_xlabel(r'$t$')
	cx2.legend(frameon=False)
	plt.show()
	'''
	
	

#Finish MPI
#MPI.COMM_WORLD.Barrier()
#MPI.Finalize()


