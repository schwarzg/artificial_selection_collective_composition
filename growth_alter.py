import numpy as np
import time
from analytic_results import *
import scipy.stats as st

def growth_sampling_m(m0,sig2m0,w0,sig2w0,covwm0,t,r,mu,s,nsample=1):
	mt=barm_th(m0,w0,t,r,mu,s)
	sig2mt=sig2m_th(m0,sig2m0,w0,sig2w0,covwm0,t,r,mu,s)
	rv=st.truncnorm(-mt/np.sqrt(sig2mt),np.inf,loc=mt,scale=np.sqrt(sig2mt))
	if nsample==1:
		return int(rv.rvs())
	else:
		return rv.rvs(size=nsmaple).astype(int)

def growth_sampling_w(w0,sig2w0,t,r,mu,s,nsample=1):
	if w0==0:
		return np.zeros(nsample)
	wt=barw_th(w0,t,r,mu)
	sig2wt=sig2w_th(w0,sig2w0,t,r,mu)
	rv=st.truncnorm(-wt/np.sqrt(sig2wt),np.inf,loc=wt,scale=np.sqrt(sig2wt))
	if nsample==1:
		return int(rv.rvs())
	else:
		return rv.rvs(size=nsmaple).astype(int)


