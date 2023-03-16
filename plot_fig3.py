import numpy as np
import matplotlib

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'DejaVu Sans:italic:bold'
plt.rcParams['mathtext.bf'] = 'DejaVu Sans:italic:bold'
plt.rcParams['axes.prop_cycle']=mpl.cycler(color=['#1f77b4', '#6600ff', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])
import itertools


from analytic_results import *
#parameter prepare
mu=1e-4
r=0.5
s=3.0e-2
N0=1000
nens=300
ncycle=1000
ncomm=10
tcycle=np.log(ncomm+1)/r

dat=np.loadtxt("data/ens/N0%s_r%s_s%s_mu%s_g%s_ncycle%d_diagram_fig2.txt"%(N0,r,s,mu,ncomm,ncycle))
with np.printoptions(precision=2):
	print(dat)

dat[dat<=0.048]=0
dat[dat>0.1]=2
dat[np.where((dat<=0.10) & (dat>0.048))]=1

print(dat)



rhats=[0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
mbars=np.array([0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000])/N0

ax=plt.axes((0.1,0.1,0.5,0.5))
frange=np.arange(0.1,1.0,0.1)
import matplotlib as mpl
oldgrey=mpl.cm.get_cmap('Greys_r')
newgrey=mpl.colors.ListedColormap(oldgrey(np.linspace(0.7,1,3)))
heatmap=ax.pcolormesh(rhats,mbars,dat,cmap=newgrey,shading='flat')
heatmap.set_clim(0,3)
cbar=plt.colorbar(heatmap,extend='max',ticks=[0,1,2])
cbar.set_label(r'Relative error $d$')
cbar.set_ticklabels([0,0.05,0.1])
ax.set_xlim(xmin=0,xmax=1)
ax.set_ylim(ymin=0,ymax=1)
ax.set_xlabel(r'Initial Frequency $\mathbf{\bar{f}_o}$',weight='bold')
ax.set_ylabel(r'Target Frequency $\mathbf{\hat{f}}$',weight='bold')

fl=0.28580445	#obtained from Fig4
fu=0.68687363
#ax.fill_between([0,fl],[fl,fl],color='lightgray',alpha=0.5)
#ax.fill_between([0,1],[1,1],[fu,fu],color='lightgray',alpha=0.5)
ax.hlines(fl,0,fl,colors='black',ls='--')
ax.vlines(fl,0,fl,colors='black',ls='--')
ax.hlines(fu,0,1,colors='black',ls='-')
ax.annotate(r'$\mathbf{f^L}$',weight='bold',xy=(0,fl),xytext=(0.10,fl+0.10),arrowprops=dict(arrowstyle='->'))
ax.annotate(r'$\mathbf{f^L}$',weight='bold',xy=(fl,0),xytext=(fl+0.10,0.10),arrowprops=dict(arrowstyle='->'))
ax.annotate(r'$\mathbf{f^U}$',weight='bold',xy=(0,fu),xytext=(0.10,fu-0.10),arrowprops=dict(arrowstyle='->'))

ax.annotate(r'Success',weight='bold',xy=(0.03,fl-0.23))
ax.annotate(r'Success',weight='bold',xy=(0.70,fu+0.10))
ax.annotate(r'Fail',weight='bold',xy=(0.45,0.35))

#Suppoting information
ax.scatter([0.05],[0.85],c='C1',marker='^')
ax.scatter([0.05],[0.50],c='C2',marker='^')
ax.scatter([0.05],[0.15],c='C3',marker='^')
ax.scatter([0.5],[0.85],c='C1',marker='v')
ax.scatter([0.5],[0.50],c='C2',marker='v')
ax.scatter([0.5],[0.15],c='C3',marker='v')


formatter='svg' #or png
filename='figures/Fig3.'
plt.savefig(filename+formatter,dpi=300,bbox_inches='tight',format=formatter)
plt.show()
