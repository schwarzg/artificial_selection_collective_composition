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

mu=1e-4
r=0.5
s=3.0e-2
N0=1000
mbar=150
ncomm=100
rhat=0.10
nsel=5
tcycle=np.log(ncomm+1)/r


#Fig2a : single trajectory
#Fig2a : extracted data of single trajectory

ncycle=13
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


cset2=[
'#1F77B4','#0D89E0','#273E4F'
]

#fl=0.28580445	#obtained from Fig4
#fu=0.68687363

ax=plt.axes((0.05,0.65,0.3,0.3)) #fhat 0.9
#ax.annotate('a',(-0.25,1.05),xycoords='axes fraction',fontweight='bold')
#bx=plt.axes((0.45,0.65,0.3,0.3)) #fhat 0.5
#bx.annotate('b',(-0.25,1.05),xycoords='axes fraction',fontweight='bold')
#cx=plt.axes((0.85,0.65,0.3,0.3)) #fhat 0.1
#cx.annotate('c',(-0.25,1.05),xycoords='axes fraction',fontweight='bold')

nens=20#300
ncycle=1000

#Read AGS data - selected only
rhat=0.1
mbars=[150]
scatterpoints=np.power(10,np.arange(0,3,0.3)).astype(int)
for idx,mbar in enumerate(mbars):

	#or, ensembled data	
	descriptor="AGS_PD_samp_N0%s_mbar%s_r%s_s%s_mu%s_ncomm%d_rhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,rhat,ncycle) 
	folder="data/ens/"+descriptor
	descriptor2="AGS_PD_samp_N0%s_mbar%s_r%s_s%s_mu%s_ncomm%d_nsel%d_rhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,nsel,rhat,ncycle) 
	folder2="data/ens/"+descriptor2

	choavg,chostd,_,_=np.loadtxt(folder+"_nens%d.cycle"%(nens),unpack=True)
	choavg2,chostd2,_,_=np.loadtxt(folder2+"_nens%d.cycle"%(nens),unpack=True)
	T=np.arange(len(choavg))+1
	
	#cx3.plot(T,choavg,marker='o',ms=3,lw=1,label=r'$\hat{f}=%s$'%rhat,c=cset2[idx])
	ax.plot(T,choavg,lw=1,c='C%d'%(idx),label='Top 1')#,c=cset2[idx])
	ax.plot(T,choavg2,lw=1,c='C%d'%(idx),ls='--',label='Top 5%')#,c=cset2[idx])
	#cx3.scatter(T[scatterpoints],choavg[scatterpoints],c='C%d'%(idx+1),marker='^',)#,c=cset2[idx])
	#cx3.fill_between(T,choavg+chostd,choavg-chostd,color='C%d'%(idx+1),alpha=0.2)
ax.hlines(rhat,0,len(T),ls=':',colors='black')#'C%d'%(4))
#ax.hlines(fu,len(T)-700,len(T),ls='-',colors='black')
#ax.hlines(fl,len(T)-700,len(T),ls='--',colors='black')
ax.set_xlim(xmin=1)
#ax.legend(fontsize='small',handlelength=1,labelspacing=0.3,loc=(0.6,0.15))
ax.set_xscale('log')
ax.legend(frameon=False,handlelength=1)
ax.set_xlabel(r'Cycle $k$')
ax.set_ylabel(r'Selected Frequency $f^*$')
ax.set_xlim(xmin=1)
'''
rhat=0.5
mbars=[0,600,950]
scatterpoints=np.power(10,np.arange(0,3,0.3)).astype(int)
for idx,mbar in enumerate(mbars):

	#or, ensembled data	
	descriptor="AGS_PD_samp_N0%s_mbar%s_r%s_s%s_mu%s_ncomm%d_rhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,rhat,ncycle) 
	folder="data/ens/"+descriptor
	choavg,chostd,_,_=np.loadtxt(folder+"_nens%d.cycle"%(nens),unpack=True)
	T=np.arange(len(choavg))+1
	
	#cx3.plot(T,choavg,marker='o',ms=3,lw=1,label=r'$\hat{f}=%s$'%rhat,c=cset2[idx])
	bx.plot(T,choavg,lw=1,label=r'$\hat{f}=%s$'%rhat,c='C%d'%(2))#,c=cset2[idx])
	#cx3.scatter(T[scatterpoints],choavg[scatterpoints],c='C%d'%(idx+1),marker='^',)#,c=cset2[idx])
	#cx3.fill_between(T,choavg+chostd,choavg-chostd,color='C%d'%(idx+1),alpha=0.2)
bx.hlines(rhat,0,len(T),ls=':',colors='C%d'%(2))
bx.hlines(fu,len(T)-700,len(T),ls='-',colors='black')
bx.hlines(fl,len(T)-700,len(T),ls='--',colors='black')

rhat=0.1
mbars=[0,200,450,950]
scatterpoints=np.power(10,np.arange(0,3,0.3)).astype(int)
for idx,mbar in enumerate(mbars):

	#or, ensembled data	
	descriptor="AGS_PD_samp_N0%s_mbar%s_r%s_s%s_mu%s_ncomm%d_rhat%s_ncycle%d"%(N0,mbar,r,s,mu,ncomm,rhat,ncycle) 
	folder="data/ens/"+descriptor
	choavg,chostd,_,_=np.loadtxt(folder+"_nens%d.cycle"%(nens),unpack=True)
	T=np.arange(len(choavg))+1
	
	#cx3.plot(T,choavg,marker='o',ms=3,lw=1,label=r'$\hat{f}=%s$'%rhat,c=cset2[idx])
	cx.plot(T,choavg,lw=1,label=r'$\hat{f}=%s$'%rhat,c='C%d'%(3))#,c=cset2[idx])
	#cx3.scatter(T[scatterpoints],choavg[scatterpoints],c='C%d'%(idx+1),marker='^',)#,c=cset2[idx])
	#cx3.fill_between(T,choavg+chostd,choavg-chostd,color='C%d'%(idx+1),alpha=0.2)
cx.hlines(rhat,0,len(T),ls=':',colors='C%d'%(3))
cx.hlines(fu,len(T)-700,len(T),ls='-',colors='black')
cx.hlines(fl,len(T)-700,len(T),ls='--',colors='black')

ax.set_ylim(ymin=0,ymax=1)
ax.set_xlabel(r'Cycle $k$')
ax.set_ylabel(r'Selected Freq. $\langle f^*\rangle$')
ax.set_xlim(xmin=1)
#ax.legend(fontsize='small',handlelength=1,labelspacing=0.3,loc=(0.6,0.15))
ax.set_xscale('log')
bx.set_xlim(xmin=1)
bx.set_xlabel(r'Cycle $k$')
bx.set_ylim(ymin=0,ymax=1)
bx.set_xlim(xmin=1)
bx.set_xscale('log')
cx.set_xlim(xmin=1)
cx.set_xlabel(r'Cycle $k$')
cx.set_ylim(ymin=0,ymax=1)
cx.set_xlim(xmin=1)
cx.set_xscale('log')

dx=plt.axes((0.10,0.0,0.5,0.5)) #fhat 0.9
dx.annotate('d',(-0.25,1.05),xycoords='axes fraction',fontweight='bold')

dat=np.loadtxt("data/ens/N0%s_r%s_s%s_mu%s_g%s_ncycle%d_diagram_fig2.txt"%(N0,r,s,mu,ncomm,ncycle))
with np.printoptions(precision=2):
	print(dat)

dat[dat<=0.048]=0
dat[dat>0.1]=2
dat[np.where((dat<=0.10) & (dat>0.048))]=1

print(dat)



rhats=[0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
mbars=np.array([0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000])/N0

#dx=plt.axes((0.1,0.1,0.5,0.5))
frange=np.arange(0.1,1.0,0.1)
#import matplotlib as mpl
oldgrey=mpl.cm.get_cmap('Greys_r')
newgrey=mpl.colors.ListedColormap(oldgrey(np.linspace(0.7,1,3)))
heatmap=dx.pcolormesh(rhats,mbars,dat,cmap=newgrey,shading='flat')
heatmap.set_clim(0,3)
cbar=plt.colorbar(heatmap,extend='max',ticks=[0,1,2])
cbar.set_label(r'Relative error $(f^*_{1000}-\hat{f})/\hat{f}$')
cbar.set_ticklabels([0,0.05,0.1])
dx.set_xlim(xmin=0,xmax=1)
dx.set_ylim(ymin=0,ymax=1)
dx.set_xlabel(r'Initial Frequency $\mathbf{\bar{f}_o}$',weight='bold')
dx.set_ylabel(r'Target Frequency $\mathbf{\hat{f}}$',weight='bold')

fl=0.28580445	#obtained from Fig4
fu=0.68687363
#dx.fill_between([0,fl],[fl,fl],color='lightgray',alpha=0.5)
#dx.fill_between([0,1],[1,1],[fu,fu],color='lightgray',alpha=0.5)
dx.hlines(fl,0,fl,colors='black',ls='--')
dx.vlines(fl,0,fl,colors='black',ls='--')
dx.hlines(fu,0,1,colors='black',ls='-')
dx.annotate(r'$\mathbf{f^L}$',weight='bold',xy=(0,fl),xytext=(0.10,fl+0.10),arrowprops=dict(arrowstyle='->'))
dx.annotate(r'$\mathbf{f^L}$',weight='bold',xy=(fl,0),xytext=(fl+0.10,0.10),arrowprops=dict(arrowstyle='->'))
dx.annotate(r'$\mathbf{f^U}$',weight='bold',xy=(0,fu),xytext=(0.10,fu-0.10),arrowprops=dict(arrowstyle='->'))

dx.annotate(r'Success',weight='bold',xy=(0.03,fl-0.23))
dx.annotate(r'Success',weight='bold',xy=(0.70,fu+0.10))
dx.annotate(r'Fail',weight='bold',xy=(0.45,0.35))
'''
'''
#Suppoting information
dx.scatter([0.05],[0.85],c='C1',marker='^')
dx.scatter([0.05],[0.50],c='C2',marker='^')
dx.scatter([0.05],[0.15],c='C3',marker='^')
dx.scatter([0.5],[0.85],c='C1',marker='v')
dx.scatter([0.5],[0.50],c='C2',marker='v')
dx.scatter([0.5],[0.15],c='C3',marker='v')
'''
'''
cxx=0.7
cxy=0.25
#Conclutsion for Two strain case
#large frequency
ex1=plt.axes((cxx,cxy,0.4,0.25))
ex1.annotate('e',(-0.05,0.8),xycoords='axes fraction',fontweight='bold')
ex1.annotate(r'Success or Fail according to the target composition',(0.00,0.8),xycoords='axes fraction')

#frame
ex1.annotate('',xy=(1,0),xytext=(0,0),arrowprops=dict(arrowstyle='->'))
ex1.fill_between([0,1],0.02,0,color='grey',alpha=0.3)
ex1.annotate(r'$\mathbf{0}$',xy=(0,-0.02))
ex1.annotate(r'$\mathbf{1}$',xy=(0.95,-0.02))

ex1.set_ylim(ymin=-0.05,ymax=0.05)


ini1=0.2
h1=0.01
ini2=0.5
h2=0.01
ini3=0.8
h3=0.01
tgt=0.9
ex1.vlines(tgt,0,0.02,colors='black',ls=':',label='Target') #target
ex1.annotate('Target',(tgt,0.02))
ex1.scatter(ini1,h1,c='black',marker='o',label='Initial') #initial 1
ex1.scatter(ini2,h2,c='black',marker='o') #initial 2
ex1.scatter(ini3,h3,c='black',marker='o') #initial 2

ex1.annotate('',xy=(ini2,h2),xytext=(ini1,h1),arrowprops=dict(arrowstyle='->'))
ex1.annotate('',xy=(ini3,h3),xytext=(ini2,h2),arrowprops=dict(arrowstyle='->'))
ex1.annotate('',xy=(tgt,h3),xytext=(ini3,h3),arrowprops=dict(arrowstyle='->'))

ex1.vlines([0.3,0.7],0.005,-0.005,colors='black')
ex1.annotate(r'$\mathbf{f^L}$',xy=(0.28,-0.02))
ex1.annotate(r'$\mathbf{f^U}$',xy=(0.68,-0.02))

ex1.set_xlim(xmin=0,xmax=1)

ex1.axis('off')

#mid frequency
ex2=plt.axes((cxx,cxy-0.15,0.4,0.25))

#frame
ex2.annotate('',xy=(1,0),xytext=(0,0),arrowprops=dict(arrowstyle='->'))
ex2.fill_between([0,1],0.02,0,color='white',alpha=0.3)
ex2.annotate(r'$\mathbf{0}$',xy=(0,-0.02))
ex2.annotate(r'$\mathbf{1}$',xy=(0.95,-0.02))

ex2.set_ylim(ymin=-0.05,ymax=0.05)

tgt=0.5
ini1=0.2
h1=0.010
ini2=0.5
h2=0.01
ini3=0.8
h3=0.01
ex2.vlines(tgt,0,0.02,colors='black',ls=':',label='Target') #target
ex2.scatter(ini1,h1,c='black',marker='o',label='Initial(low)') #initial 1
ex2.scatter(ini2,h2,c='black',marker='o',label='Initial(mid)') #initial 2
ex2.scatter(ini3,h3,c='black',marker='o',label='Initial(high)') #initial 2

ex2.annotate('',xy=(ini2,h2),xytext=(ini1,h1),arrowprops=dict(arrowstyle='->'))
ex2.annotate('',xy=(0.7,h2),xytext=(ini2,h2),arrowprops=dict(arrowstyle='->'))
ex2.annotate('',xy=(0.7,h3),xytext=(ini3,h3),arrowprops=dict(arrowstyle='->'))

ex2.vlines([0.3,0.7],0.005,-0.005,colors='black')
ex2.annotate(r'$\mathbf{f^L}$',xy=(0.28,-0.02))
ex2.annotate(r'$\mathbf{f^U}$',xy=(0.68,-0.02))

ex2.set_xlim(xmin=0,xmax=1)
ex2.axis('off')

#small frequency
ex3=plt.axes((cxx,cxy-0.3,0.4,0.25))

#frame
ex3.annotate('',xy=(1,0),xytext=(0,0),arrowprops=dict(arrowstyle='->'))
ex3.fill_between([0,0.3],0.025,0,color='grey',alpha=0.3, label='success')
ex3.fill_between([0.3,1],0.025,0,color='white',alpha=0.3,label='fail')
ex3.annotate(r'$\mathbf{0}$',xy=(0,-0.02))
ex3.annotate(r'$\mathbf{1}$',xy=(0.95,-0.02))
ex3.annotate('Mutant frequency',xy=(0.5,0.2),xycoords='axes fraction',ha='center',va='top')

ex3.set_ylim(ymin=-0.05,ymax=0.05)

tgt=0.05
ini1=0.2
h1=0.015
ini2=0.5
h2=0.015
ini3=0.8
h3=0.015

ex3.scatter(ini1,h1,c='black',marker='o',label='Initial') #initial 1
ex3.scatter(ini2,h2,c='black',marker='o') #initial 2
ex3.scatter(ini3,h3,c='black',marker='o') #initial 2
ex3.vlines(tgt,-0.005,0.025,colors='black',ls=':',label='Target') #target

ex3.annotate('',xy=(tgt,h1),xytext=(ini1,h1),arrowprops=dict(arrowstyle='->'))
ex3.annotate('',xy=(0.7,h2),xytext=(ini2,h2),arrowprops=dict(arrowstyle='->'))
ex3.annotate('',xy=(0.7,h3),xytext=(ini3,h3),arrowprops=dict(arrowstyle='->'))

ex3.vlines([0.3,0.7],0.005,-0.005,colors='black')
ex3.annotate(r'$\mathbf{f^L}$',xy=(0.28,-0.02))
ex3.annotate(r'$\mathbf{f^U}$',xy=(0.68,-0.02))

#ex3.legend(frameon=False,loc=(1.1,0.5))
ex3.set_xlim(xmin=0,xmax=1)
ex3.axis('off')
'''
formatter='svg' #or 'png'
#plt.savefig('figures/FigS7.'+formatter,dpi=300,bbox_inches='tight',format=formatter)
plt.show()

