import numpy as np
import matplotlib.pyplot as plt

#parameter prepare
mu=1e-4
r=0.5
s=3e-2
N0=900
mbar=100
ncomm=10
rhat=0 if s>0 else 1
nens=300
ncycle=10
tcycle=np.log(ncomm+1)/r

folder="data/cond/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(N0,r,mu,s,ncomm,nens)
if not os.path.exists(folder):
	os.mkdir(folder)

import scipy.stats as st
import scipy.integrate as intg
from scipy.integrate import simpson
import scipy.special as spc
from analytic_results import *

from one_step_functions import *

from model import model
from tau_leaping import tau_leaping
from direct import direct
from growth_alter import *
from selection import select
from selection import community_function as cf
from selection import score_function as sf
from reproduction import hypergeometric_reproduce
from stats_custom_tools import * 

#f0=0.1
#frequency loop
fdat=[]
fseldat=[]
fextdat=[]
means=[]
medians=[]
mediansex=[]
meansex=[]
for f0 in np.arange(0.00,0.99,0.01): 
#for f0 in np.arange(0.00,0.99,0.1): 
	print("Frequency ",f0)
	fbins=np.linspace(np.maximum(0,f0-0.05),np.minimum(1,f0+0.05),30)
	#m0=int(f0*Nf)
	#w0=Nf-m0
	#print(w0,m0,f0)
	############################
	# Simulation part
	############################
	#growth phase
	#prepare data space
	#wsel=np.array([])
	#msel=np.array([])
	w_sel=np.zeros(nens)
	m_sel=np.zeros(nens)
	fsampmean=np.array([])
	fsampstd=np.array([])
	
	#ensmble loop
	dats=[]
	for e in range(nens):
		#print("ensemble ",e)
		#reproduction phase	
		m0sel=np.random.binomial(N0,f0,size=ncomm)	
		m0sel=np.sort(m0sel)
		w0sel=N0-m0sel
		
		#For each community
		lastw=np.zeros(ncomm)
		lastm=np.zeros(ncomm)
		for j in range(ncomm):
			#Growth phase
			wsamp=growth_sampling_w(w0sel[j],0,tcycle,r,mu,s)
			msamp=growth_sampling_m(m0sel[j],0,w0sel[j],0,0,tcycle,r,mu,s)
	
			#store for selection
			lastw[j]=wsamp
			lastm[j]=msamp
		
		#Selection phase
		ind_sel=select(lastw,lastm,r=r,s=s,rhat=rhat)
		w_sel[e]=lastw[ind_sel]
		m_sel[e]=lastm[ind_sel]
		
	cfs_sel=cf(w_sel,m_sel)
	dat,_,_=plt.hist(cfs_sel,bins=fbins,density=True)
	fseldat.append(dat)	#store histogram of selected data
	medians.append(np.median(cfs_sel))
	means.append(np.mean(cfs_sel))
	
	############################
	# Evaluation part
	############################	
	pdf=transition_probability_iid_pdf(fbins[:-1],f0,tcycle,N0,r,mu,s)
	pdfnorm=simpson(pdf,dx=fbins[1]-fbins[0])
	pdf=pdf/pdfnorm
	cdf=np.cumsum(pdf)*(fbins[1]-fbins[0])
	cdf=cdf/cdf[-1]
	extpdf=ncomm*(0.5+np.sign(s)*(0.5-cdf))**(ncomm-1)*pdf #Same dimension with dat
	#print(pdf,pdfnorm)
	#print(cdf)
	#plt.plot(fbins[:-1],pdf)
	#plt.plot(fbins[:-1],cdf)
	#plt.ylim(ymin=0,ymax=1)
	fextdat.append(extpdf)	
	median=quantile(0.5,fbins,extpdf,x0=np.min(fbins))
	mediansex.append(median)
	mean=np.sum(fbins[:-1]*extpdf*(fbins[1]-fbins[0]))
	meansex.append(mean)
#np.savetxt("data/one/one_step_transition_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s_fine.dat_v2"%(N0,r,mu,s,ncomm,nens),fdat)
np.savetxt(folder+"_sel.dat",fseldat)
np.savetxt(folder+"_smedian.dat",medians)
np.savetxt(folder+"_smean.dat",means)
np.savetxt(folder+"_ext.dat",fextdat)
np.savetxt(folder+"_emedian.dat",mediansex)
np.savetxt(folder+"_emean.dat",meansex)
print(medians)
plt.show()	
