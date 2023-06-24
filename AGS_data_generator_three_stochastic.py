'''

Additional model for two mutation system

	w -r-> w+w
	m -r+s-> m+m
	v -r+2s-> v+v
	w -mu-> m -mu-> v

	target composition

	1. fm only (line target)
	2. fm and fv (point target)

'''

import numpy as np
from gillespie.model import model
from gillespie.tau_leaping import tau_leaping
import os,sys
from tqdm import tqdm
#import matplotlib.pyplot as plt

def cf1(w,m,v):
	#community function for m 
	N=w+m+v
	f=np.divide(m,N)
	return f

def cf2(w,m,v):
	#community function for v
	N=w+m+v
	f=np.divide(v,N)
	return f

def select1(w,m,v,mhat):
	comms=cf1(w,m,v)
	d=np.fabs(mhat-comms)
	smindist=np.argmin(d)
	return smindist
def select2(w,m,v,mhat):
	f1=cf1(w,m,v)
	f2=cf2(w,m,v)
	d=(mhat-f1)**2+(vhat-f2)**2
	smindist=np.argmin(d)
	return smindist

def reproduce(w0,m0,v0,N0,Ngroup=None):
	N=w0+m0+v0
	if Ngroup==None:
		Ng=int(N/N0)
	else:
		Ng=Ngroup

	Wnew=[]
	Mnew=[]	
	Vnew=[]
	N=w0+m0+v0
	fm,fv=cf1(w0,m0,v0),cf2(w0,m0,v0)
	for j in range(Ng):
		if w0>=0 and m0>=0:
			if N0<w0+m0:
				msamp,vsamp,wsamp=np.random.multinomial(N0,[fm,fv])	
				Wnew.append(wsamp)
				Mnew.append(msamp)
				Vnew.append(vsamp)
			else:
				Mnew.append(m0)
				Wnew.append(w0)
				Vnew.append(v0)
			m0=m0-Mnew[-1]
			w0=w0-Wnew[-1]
			v0=w0-Vnew[-1]
		else:
			print("Somepart is negative")
			return -1,-1
	return np.array(Wnew),np.array(Mnew),np.array(Vnew)

mu=1e-4
r=0.5
s=3e-2
N0=1000
m0=85
v0=5
ncomm=10
mhat=0.02
vhat=0.9
nensemble=30
ncycle=30
tcycle=np.log(ncomm+1)/r

if len(sys.argv)!=1:
	m0=int(sys.argv[1])	
	v0=int(sys.argv[2])	
	mhat=float(sys.argv[3])
	vhat=float(sys.argv[4])

descriptor="AGS_PD_sto_N0%s_mbar%s_vbar%s_r%s_s%s_mu%s_ncomm%d_mhat%s_vhat%s_ncycle%d"%(N0,m0,v0,r,s,mu,ncomm,mhat,vhat,ncycle)
folder="data/raw/"+descriptor

'''
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')
parser.add_argument('--sum', dest='accumulate', action='store_const',
                    const=sum, default=max,
                    help='sum the integers (default: find the max)')
'''

proFunc=[
lambda x: r*x[0],
lambda x: mu*x[0],
lambda x: (r+s)*x[1],
lambda x: mu*x[1],
lambda x: (r+2*s)*x[2]
]
changeVec=np.array([[1,0,0],[-1,1,0],[0,1,0],[0,-1,1],[0,0,1]]).astype(int)
rxnOrder=np.array([[1,0,0],[1,0,0],[0,1,0],[0,1,0],[0,0,1]]).astype(int)

growthmodel=model(proFunc=proFunc,changeVec=changeVec,rxnOrder=rxnOrder)

solver=tau_leaping(model=growthmodel,eps=0.01)

M0=np.ones(ncomm,dtype=int)*m0
V0=np.ones(ncomm,dtype=int)*v0
W0=N0-M0-V0

nproc=1
rank=0

for e in range(rank,nensemble,nproc):
#for e in tqdm(range(rank,nensemble,nproc)):
	#print('ens ',e)	
	W0sel=np.copy(W0)
	M0sel=np.copy(M0)
	V0sel=np.copy(V0)
	
	#prepare data space
	tsel=np.array([])
	wsel=np.array([])
	msel=np.array([])
	vsel=np.array([])

		
	t_selected=np.array([])
	w_selected=np.array([])
	m_selected=np.array([])
	v_selected=np.array([])
	ind_selected=np.array([])
	
	#run community lifecycle
	for i in tqdm(range(ncycle),desc='ens %d'%e):
		#print('cycle',i)
		w_selected=np.append(w_selected,W0sel)
		m_selected=np.append(m_selected,M0sel)
		v_selected=np.append(v_selected,V0sel)
		
		#Space for selection	
		lastw=np.array([])
		lastm=np.array([])
		lastv=np.array([])

		#For each community
		for j in range(ncomm):

			#Growth phase
			T,X=solver.run(np.array([W0sel[j],M0sel[j],V0sel[j]]),0,tcycle)
			tsel=np.append(tsel,T[:-1]+i*tcycle)	
			wsel=np.append(wsel,X[0,:-1])
			msel=np.append(msel,X[1,:-1])
			vsel=np.append(vsel,X[2,:-1])

			#store for selection
			lastw=np.append(lastw,X[0,-1])
			lastm=np.append(lastm,X[1,-1])	
			lastv=np.append(lastv,X[2,-1])	
			'''	
			#Trajectory data temporary test
			if e==0:
				import matplotlib.pyplot as plt
				#cfs=cf(X[0],X[1])
				#plt.plot(T+i*tcycle,cfs,c=cset[j])
			'''

		#Selection phase
		w_selected=np.append(w_selected,lastw)
		m_selected=np.append(m_selected,lastm)
		v_selected=np.append(v_selected,lastv)
		#ind_sel=select1(lastw,lastm,lastv,mhat)
		ind_sel=select2(lastw,lastm,lastv,mhat)
		w_sel=lastw[ind_sel]
		m_sel=lastm[ind_sel]
		v_sel=lastv[ind_sel]
		ind_selected=np.append(ind_selected,ind_sel)
		t_selected=np.append(t_selected,(i+1)*tcycle)	
		
		#reproduction phase	
		#w0sel,m0sel,v0sel=reproduce(w_sel,m_sel,v_sel,N0,Ngroup=10)
		fm=cf1(w_sel,m_sel,v_sel)
		fv=cf2(w_sel,m_sel,v_sel)
		#print(e,i,fm,fv,fm+fv,1-fm-fv<0)
		if 1-fm-fv<0:
			X=np.zeros((ncomm,3))
			#X[:,0]=np.zeros(ncomm)
			X[:,1]=np.random.binomial(N0,fm,size=ncomm)	
			X[:,2]=N0-X[:,1]
			print('No w')	
		else:
			X=np.random.multinomial(N0,[1-fm-fv,fm,fv],size=ncomm)	
		W0sel=X[:,0]
		M0sel=X[:,1]
		V0sel=X[:,2]

	#Trajectory data temporary test
	#if e==0:
		#cfs=cf(X[0],X[1])
		#plt.plot(T+i*tcycle,cfs,c=cset[j])
		#plt.plot(np.arange(ncycle),w_selected[::2*ncomm],c='C0')
		#plt.plot(np.arange(ncycle),m_selected[::2*ncomm],c='C1')
		#plt.plot(np.arange(ncycle),v_selected[::2*ncomm],c='C2')
			
	#Save data
	#np.savetxt for data at cycle time only
	output=np.hstack((np.array([t_selected,ind_selected]).T,w_selected.reshape(ncycle,2*ncomm),m_selected.reshape(ncycle,2*ncomm),v_selected.reshape(ncycle,2*ncomm)))
	np.savetxt(folder+"%d.cycle"%(e),output)

#plt.show()
