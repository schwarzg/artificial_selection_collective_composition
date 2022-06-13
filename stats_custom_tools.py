'''

	My custom functions for statistics

'''

import numpy as np
import scipy.stats as st
import scipy.interpolate as itp
import scipy.optimize as opt
import matplotlib.pyplot as plt

#quartile of 1D histogram
def quantile(q,bins,hist,x0=0,qmin=0,qmax=1):
	x=bins[:-1]+0.5*(bins[1]-bins[0])
	cums=np.cumsum(hist*(bins[1]-bins[0]))
	y=itp.interp1d(x,cums-q,fill_value='extrapolate')
	xp=np.linspace(bins[0],x[-1],100)
	#sol=opt.root(y,x0,method='broyden2').x
	sol=opt.brentq(y,x[0],x[-1])
	return sol
	
	
