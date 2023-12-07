import numpy as np
from model import model
from tau_leaping import tau_leaping
from selection import select
from selection import community_function as cf
from selection import score_function as sf
import matplotlib as mpl
import matplotlib.pyplot as plt
import itertools as itt

plt.rcParams["font.family"] = "arial"
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
#plt.rcParams["mathtext.fontset"]='stix'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'STIXGeneral:italic:bold'
plt.rcParams['mathtext.bf'] = 'STIXGeneral:italic:bold'
plt.rcParams['axes.prop_cycle']=mpl.cycler(color=['#1f77b4', '#6600ff', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])
plt.rcParams['axes.prop_cycle']=mpl.cycler(color=['#1f77b4', '#6600ff', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])

mu=1e-4
r=0.5
s=3.0e-2
N0=1000
mbar=200
ncomm=10
rhat=0.10

tcycle=np.log(ncomm+1)/r


#Fig2a : single trajectory
#Fig2a : extracted data of single trajectory

ncycle=13
cx=plt.axes((0.10,0.55,0.33,0.33))
################################
#data generate version
################################

m0sel=np.sort(np.clip(np.random.binomial(N0,mbar/N0,size=ncomm).astype(int),0,1000))
w0sel=N0-m0sel

m0nat=m0sel
w0nat=w0sel

#define model
proFunc=[
lambda x: r*x[0],
lambda x: mu*x[0],
lambda x: (r+s)*x[1]
]
changeVec=np.array([[1,0],[-1,1],[0,1]]).astype(int)
rxnOrder=np.array([[1,0],[1,0],[0,1]]).astype(int)
growthmodel=model(proFunc=proFunc,changeVec=changeVec,rxnOrder=rxnOrder)
solver=tau_leaping(model=growthmodel,eps=0.01)

selected=[]
averaged=np.array([])
averagedGR=np.array([])
deviated=np.array([])
deviatedGR=np.array([])
minGR=np.array([])
maxGR=np.array([])
tselect=[]

for i in range(ncycle):
	print("Cycle %d"%i)
	tsel=[] #group,time
	wsel=[]
	msel=[]
	selind=[]
	tnat=[]
	wnat=[]
	mnat=[]
	for j in range(ncomm):
		tsel.append(np.array([]))
		wsel.append(np.array([]))
		msel.append(np.array([]))
		tnat.append(np.array([]))
		wnat.append(np.array([]))
		mnat.append(np.array([]))

	#growth phase	
	lastw=[]
	lastm=[]
	lastwGR=[]
	lastmGR=[]
	for j in range(ncomm):

		#ACS model : Growth, Selection, and Reproduction
		T,X=solver.run(np.array([w0sel[j],m0sel[j]]),0,tcycle)
		tsel[j]=np.append(tsel[j],T[:-1]+i*tcycle)	
		wsel[j]=np.append(wsel[j],X[0,:-1])
		msel[j]=np.append(msel[j],X[1,:-1])
	
		#store for selection
		lastw.append(X[0,-1])
		lastm.append(X[1,-1])

		#NS: Growth and Reproduction without selection
		T,X=solver.run(np.array([w0nat[j],m0nat[j]]),0,tcycle)
		tnat[j]=np.append(tnat[j],T[:-1]+i*tcycle)	
		wnat[j]=np.append(wnat[j],X[0,:-1])
		mnat[j]=np.append(mnat[j],X[1,:-1])
		
		lastwGR.append(X[0,-1])
		lastmGR.append(X[1,-1])
		
		m0nat[j]=np.random.binomial(N0,cf(wnat[j][-1],mnat[j][-1],r=r,s=s),size=1)
		w0nat[j]=N0-m0nat[j]
	
	lastw=np.array(lastw)
	lastm=np.array(lastm)
	lastwGR=np.array(lastwGR)
	lastmGR=np.array(lastmGR)
	cfs=cf(lastwGR,lastmGR)
	averagedGR=np.append(averagedGR,np.mean(cfs))
	deviatedGR=np.append(deviatedGR,np.std(cfs))
	maxGR=np.append(maxGR,np.max(cfs))
	minGR=np.append(minGR,np.min(cfs))
	
	#selection phase
	ind_sel=select(lastw,lastm,r=r,s=s,rhat=rhat)
	w_sel=lastw[ind_sel]
	m_sel=lastm[ind_sel]
	selind.append(ind_sel)
	selected.append(cf(w_sel,m_sel,r=r,s=s))
	cfs=cf(lastw,lastm)
	averaged=np.append(averaged,np.mean(cfs))
	deviated=np.append(deviated,np.std(cfs))
	tselect.append((i+1)*tcycle)

	#save data	
	fts=open('data/figS/fig2_tsel_%d.dat'%i,'a')
	fws=open('data/figS/fig2_wsel_%d.dat'%i,'a')
	fms=open('data/figS/fig2_msel_%d.dat'%i,'a')
	fsi=open('data/figS/fig2_selind_%d.dat'%i,'a')
	np.savetxt(fsi,selind)
	ftn=open('data/figS/fig2_tnat_%d.dat'%i,'a')
	fwn=open('data/figS/fig2_wnat_%d.dat'%i,'a')
	fmn=open('data/figS/fig2_mnat_%d.dat'%i,'a')
	for j in range(ncomm):
		np.savetxt(fts,tsel[j],newline=' ')
		fts.write('\n')
		np.savetxt(fws,wsel[j],newline=' ')
		fws.write('\n')
		np.savetxt(fms,msel[j],newline=' ')
		fms.write('\n')
		np.savetxt(ftn,tnat[j],newline=' ')
		ftn.write('\n')
		np.savetxt(fwn,wnat[j],newline=' ')
		fwn.write('\n')
		np.savetxt(fmn,mnat[j],newline=' ')
		fmn.write('\n')
	fts.close()
	fws.close()
	fms.close()
	fsi.close()
	ftn.close()
	fwn.close()
	fmn.close()
	
	#draw	
	for j in range(0,ncomm):
		cfs=cf(wsel[j],msel[j],r=r,s=s)
		if j==ind_sel:
			cx.plot(tsel[j],cfs,lw=1,c='black')
		else:
			cx.plot(tsel[j],cfs,lw=1,c='black',alpha=0.3)
		#nx.plot(tsel[j],(wsel[j]+msel[j]),lw=1,c='red')
	
	#reproduction phase	
	#w0sel,m0sel=hypergeometric_reproduce(w_sel,m_sel,N0)	
	m0sel=np.random.binomial(N0,cf(w_sel,m_sel,r=r,s=s),size=ncomm)
	m0sel=np.sort(m0sel)
	w0sel=N0-m0sel


'''
################################
#data read version
################################
tsels=[] #group,time
wsels=[]
msels=[]
selinds=[]
tnats=[]
wnats=[]
mnats=[]

cx=plt.axes((0.10,0.55,0.33,0.33))
cx.annotate('a',(-0.25,1.05),xycoords='axes fraction',fontweight='bold')
selected=[]
averaged=np.array([])
averagedGR=np.array([])
for i in range(ncycle):
	print("Cycle %d"%i)
	tsel=[] #group,time
	wsel=[]
	msel=[]
	selind=[]
	tnat=[]
	wnat=[]
	mnat=[]

	fts=open('data/figS/fig2_tsel_%d.dat'%i,'r')
	lws=open('data/figS/fig2_wsel_%d.dat'%i,'r')
	fms=open('data/figS/fig2_msel_%d.dat'%i,'r')
	fsi=open('data/figS/fig2_selind_%d.dat'%i,'r')
	ind_sel=np.loadtxt(fsi).astype(int)
	ftn=open('data/figS/fig2_tnat_%d.dat'%i,'r')
	fwn=open('data/figS/fig2_wnat_%d.dat'%i,'r')
	fmn=open('data/figS/fig2_mnat_%d.dat'%i,'r')

	#read 
	for j in range(ncomm):
		tsel.append(np.array(fts.readline().split(' ')[:-1]).astype(float))
		wsel.append(np.array(fws.readline().split(' ')[:-1]).astype(float))
		msel.append(np.array(fms.readline().split(' ')[:-1]).astype(float))
		tnat.append(np.array(ftn.readline().split(' ')[:-1]).astype(float))
		wnat.append(np.array(fwn.readline().split(' ')[:-1]).astype(float))
		mnat.append(np.array(fmn.readline().split(' ')[:-1]).astype(float))
	
	#prepare cx2
	selected.append(cf(wsel[ind_sel][-1],msel[ind_sel][-1],r=r,s=s))
	wnatnp=np.array([])
	mnatnp=np.array([])
	for j in range(ncomm):
		wnatnp=np.append(wnatnp,wnat[j][-1])
		mnatnp=np.append(mnatnp,mnat[j][-1])
	cfs=cf(wnatnp,mnatnp,r=r,s=s)
	averagedGR=np.append(averagedGR,np.mean(cfs))

	#draw	
	for j in range(0,ncomm):
		cfs=cf(wsel[j],msel[j],r=r,s=s)
		if j==ind_sel:
			cx.plot(tsel[j],cfs,lw=1,c='black')
		else:
			cx.plot(tsel[j],cfs,lw=0.5,c='gray',alpha=0.4)
'''

cx.annotate('a',(-0.25,1.05),xycoords='axes fraction',fontweight='bold')
cx.hlines(rhat,0,ncycle*tcycle,linestyles=':',colors='black')
cx.set_ylabel(r'F frequency $f$')
cx.set_xlabel(r'Time $t$')
cx.set_ylim(ymin=0.02,ymax=0.25)
cx.annotate(r'$\tau\approx4.8$',xy=(8*tcycle,0.06))
cx.annotate('',xy=(9*tcycle,0.10),xytext=(10*tcycle,0.10),arrowprops=dict(arrowstyle='->',mutation_scale=10))
cx.annotate('',xy=(8*tcycle,0.10),xytext=(7*tcycle-0.01,0.10),arrowprops=dict(arrowstyle='->',mutation_scale=10))
cx.vlines([8*tcycle,9*tcycle],0.09,0.11,linestyles='-',colors='black')

cx2=plt.axes((0.55,0.55,0.33,0.33))
cx2.annotate('b',(-0.25,1.05),xycoords='axes fraction',fontweight='bold')
cx2.set_ylabel(r'F frequency $f$')
cx2.set_xlabel(r'Cycle $k$')
cx2.hlines(rhat,0,ncycle,linestyles=':',colors='black')
k=np.arange(1,ncycle+1)
cx2.plot(k,selected,marker='^',ms=4,c='black',label=r'$f^*$')
cx2.plot(k,averagedGR,marker='s',ms=4,c='blue',ls='--',label=r'$\bar{f}_\mathrm{NS}$')
cx2.legend(fontsize='small',handlelength=1,labelspacing=0.3)
cx2.set_ylim(ymin=0.05)

formatter='svg' #or 'png'
#plt.savefig('figures/FigS4.'+formatter,dpi=300,bbox_inches='tight',format=formatter)
plt.show()

