import numpy as np
import matplotlib.pyplot as plt
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

rhats=[0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
mbars=np.array([0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000])/N0

ax=plt.axes((0.1,0.1,0.5,0.5))
frange=np.arange(0.1,1.0,0.1)
heatmap=ax.pcolormesh(rhats,mbars,dat,cmap='Blues_r',shading='flat')
heatmap.set_clim(0,0.05)
cbar=plt.colorbar(heatmap,extend='max')
cbar.set_label(r'Relative error $d$')
ax.set_xlim(xmin=0,xmax=1)
ax.set_ylim(ymin=0,ymax=1)
ax.set_xlabel(r'Initial Frequency $\bar{f}_o$')
ax.set_ylabel(r'Target Frequency $\hat{f}$')

#Suppoting information
ax.scatter([0.05],[0.85],c='C1',marker='^')
ax.scatter([0.05],[0.50],c='C2',marker='^')
ax.scatter([0.05],[0.15],c='C3',marker='^')
ax.scatter([0.5],[0.85],c='C1',marker='v')
ax.scatter([0.5],[0.50],c='C2',marker='v')
ax.scatter([0.5],[0.15],c='C3',marker='v')

#computed from Fig 4 in Main text
fl=0.28580445
fu=0.68687363
ax.hlines(fl,0,fl,colors='black',ls='--')
ax.vlines(fl,0,fl,colors='black',ls='--')
ax.hlines(fu,0,1,colors='black',ls='--')
ax.annotate(r'$f^L$',xy=(0,fl),xytext=(0.10,fl+0.10),arrowprops=dict(arrowstyle='->'))
ax.annotate(r'$f^L$',xy=(fl,0),xytext=(fl+0.10,0.10),arrowprops=dict(arrowstyle='->'))
ax.annotate(r'$f^U$',xy=(0,fu),xytext=(0.10,fu-0.10),arrowprops=dict(arrowstyle='->'))

filename='figures/FigS3.'
ftype='svg'
plt.savefig(filename+ftype,dpi=300,bbox_inches='tight',format=ftype)
#plt.savefig("Fig3_s%s.png"%s,dpi=300,bbox_inches='tight')#,format='svg')
plt.show()
