import numpy as np
import matplotlib.pyplot as plt

#parameter prepare
mu=1e-4
r=0.5
s=3e-2
N0=1000
mbar=100
ncomm=10
rhat=0 if s>0 else 1
nens=300
tcycle=np.log(ncomm+1)/r

folder="data/cond/conditional_probability_N0%s_f0_r%s_mu%s_s%s_g%s_nens%s"%(N0,r,mu,s,ncomm,nens)

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
from stats_custom_tools import * 

fdat=[]
fseldat=[]
fextdat=[]
means=[]
medians=[]
mediansex=[]
meansex=[]
f0s=np.arange(0.00,0.99,0.01)
for f0 in f0s:
	print("Frequency ",f0)
	fbins=np.linspace(np.maximum(0,f0-0.05),np.minimum(1,f0+0.05),30)
	'''	
	############################
	# Simulation part
	############################
	#prepare data space
	w_sel=np.zeros(nens)
	m_sel=np.zeros(nens)
	fsampmean=np.array([])
	fsampstd=np.array([])
	
	#ensmble loop
	dats=[]
	for e in range(nens):
		
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
	'''
	############################
	# Evaluation part
	############################	
	pdf=transition_probability_iid_pdf(fbins[:-1],f0,tcycle,N0,r,mu,s)
	pdfnorm=simpson(pdf,dx=fbins[1]-fbins[0])
	pdf=pdf/pdfnorm
	cdf=np.cumsum(pdf)*(fbins[1]-fbins[0])
	cdf=cdf/cdf[-1]
	extpdf=ncomm*(0.5+np.sign(s)*(0.5-cdf))**(ncomm-1)*pdf #Same dimension with dat
	fextdat.append(extpdf)	
	median=quantile(0.5,fbins,extpdf,x0=np.min(fbins))
	mediansex.append(median)
	mean=np.sum(fbins[:-1]*extpdf*(fbins[1]-fbins[0]))
	meansex.append(mean)
print(np.array(mediansex)-f0s)
#np.savetxt(folder+"_sel.dat",fseldat)
#np.savetxt(folder+"_smedian.dat",medians)
#np.savetxt(folder+"_smean.dat",means)
np.savetxt(folder+"_ext.dat",fextdat)
np.savetxt(folder+"_emedian.dat",mediansex)
np.savetxt(folder+"_emean.dat",meansex)
plt.show()	
