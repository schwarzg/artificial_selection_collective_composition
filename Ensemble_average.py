import numpy as np
import matplotlib.pyplot as plt
from selection import community_function as cf
import itertools as itt
import tqdm.contrib.itertools as titt
	
mu=0#1e-4
r=0.5
s=3.0e-2
N0=1000
ncycle=1000
ncomm=10

#AGS plot
dats=np.array([])	#data for saving 
rhats=[0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]
mbars=[0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950]

#rhats=[0.1]
#mbars=[150,800]

dats=[]

nensemble=100#20#300
c=0
#for (rhat,mbar) in itt.product(rhats,mbars):
for (rhat,mbar) in titt.product(rhats,mbars):
	print(rhat,mbar)	
	descriptor="AGS_PD_samp_N0%s_mbar%s_r%s_s%s_mu%s_ncomm%d_rhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,rhat,ncycle) 
	infolder="data/raw/"+descriptor+"/"
	outfolder="data/ens/"+descriptor
	AGS=[]
	for e in range(nensemble):
		AGS.append(np.loadtxt(infolder+"%d.cycle"%(e)))
	AGS=np.array(AGS) # [ensemble, timestamp, (T,selected index,w----,m----)]
	
	T=AGS[0,:,0]
	sel_inds=AGS[:,:,1].astype(int) #selected index - [ensemble, timestamp]
	c_sel=np.zeros((nensemble,len(T),ncomm)) #c avg datas, [ensemble, timestamp]
	w_cho=np.zeros((nensemble,len(T))) #w datas, [ensemble, timestamp]
	m_cho=np.zeros((nensemble,len(T))) #m datas, [ensemble, timestamp]
	#ensemble loop
	for i in range(nensemble):
		#timestamp loop
		for j in range(len(T)):
			c_sel[i,j,:]=cf(AGS[i,j,ncomm+2:ncomm*2+2],AGS[i,j,ncomm*3+2:])
			w_cho[i,j]=AGS[i,j,sel_inds[i,j]+2+ncomm]
			m_cho[i,j]=AGS[i,j,sel_inds[i,j]+2+ncomm*3]
	
	#Artificial group selection with the selected
	c_cho=cf(w_cho,m_cho)
	choavg=np.mean(c_cho,axis=0)	
	chostd=np.std(c_cho,axis=0)	
	cho75=np.quantile(c_cho,0.75,axis=0)
	cho25=np.quantile(c_cho,0.25,axis=0)
	np.savetxt(outfolder+"_nens%d.cycle"%nensemble,np.array([choavg,chostd,cho75,cho25]).T)
	'''
	choavg,chostd,cho75,cho25=np.loadtxt(outfolder+"_nens%d.cycle"%nensemble,unpack=True)
	'''	
	#save data for fig3
	dats.append(np.fabs(np.mean(choavg[-1])-rhat))#/(rhat+1e-12))

dats=np.array(dats).reshape((len(rhats),len(mbars)))
print(dats)
np.savetxt("data/ens/N0%s_r%s_s%s_mu%s_g%s_ncycle%d_diagram_abs.txt"%(N0,r,s,mu,ncomm,ncycle),dats)
