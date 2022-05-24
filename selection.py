import numpy as np


def community_function(w,m,r=0.5,s=5e-2):
	N=w+m
	f=np.divide(m,N)
	#return r+s*f
	return f

def score_function(r,rhat):
	return 1/(np.fabs(r-rhat)+1e-16)
	#return np.exp(-np.fabs(r-rhat))

def select(w,m,r=0.5,s=5e-2,rhat=0.05):
	comms=community_function(w,m,r=r,s=s)
	score=score_function(comms,rhat)
	smax=np.where(score==np.max(score))[0][0]
	#print("function:",score,np.max(score),smax)
	#return w[smax],m[smax]
	return smax

if __name__=='__main__':
	mu=1e-4
	r=0.5
	s=5e-2
	w0=1000
	m0=0
	tmax=10
	nsample=100000
	ngroup=2000
	npergroup=int(nsample/ngroup)
	folder="data/time/%d_%d_%s_%s_%s_%s/"%(w0,m0,r,s,mu,tmax)
	t=9.9
	w,m=np.loadtxt(folder+"/%.1f_tau.dat"%t,unpack=True)
	rhat=0.05
	sc=score_function(community_function(w,m),rhat)
	print("w,m,community",w[:10],m[:10],community_function(w[:10],m[:10]))
	print(select(w[:10],m[:10])	)
	
	import matplotlib.pyplot as plt
	#plt.hist(sc,bins=100,density=True)
	#plt.show()
	
	s_max=np.array([])
	for i in range(ngroup):
		s_max=np.append(s_max,np.max(sc[i*npergroup:(i+1)*npergroup]))

	import scipy.stats as st

	loc1,scale1=st.gumbel_r.fit(s_max)	
	c2,loc2,scale2=st.invweibull.fit(s_max)	
	#c3,loc3,scale3=st.weibull_min.fit(s_max)	
	c4,loc4,scale4=st.genextreme.fit(s_max)	
	
	bins=np.arange(np.min(s_max)*0.8,np.max(s_max)*1.2,(np.max(s_max)-np.min(s_max))/40)
	ax=plt.axes((0.1,0.1,0.50,0.4))
	ax.hist(s_max,density=True,bins=bins,histtype=u'step',label=r'Data,$%.2f,%.2f$'%(np.mean(s_max),np.std(s_max)))
	ax.plot(bins,st.gumbel_r.pdf(bins,loc=loc1,scale=scale1),label=r'Gumbel,$%.2f,%.2f$'%(loc1,scale1))
	#ax.plot(bins,st.invweibull.pdf(bins,c=c2,loc=loc2,scale=scale2),label=r'Frechet,$%.2f,%.2f,%.2f$'%(c2,loc2,scale2))
	#ax.plot(bins,st.weibull_min.pdf(bins,c=c3,loc=loc3,scale=scale3),label=r'Weibull,$%.2f,%.2f,%.2f$'%(c3,loc3,scale3))
	ax.plot(bins,st.genextreme.pdf(bins,c=c4,loc=loc4,scale=scale4),label=r'GEV,$%.2f,%.2f,%.2f$'%(c4,loc4,scale4))
	ax.hlines(1/100000,bins[0],bins[-1],colors='black',ls=':')
	ax.legend(frameon=False)
	ax.set_yscale('log')
	ax.set_ylim(ymin=1/100000)#,ymax=10)	
	dx=plt.axes((0.1,0.6,0.50,0.4))
	dx.hist(score_function(community_function(w,m),rhat),bins=30,density=True)
	dx.set_xlabel(r'$s$')
	dx.set_ylabel(r'$P(s)$')
	#plt.savefig("%d_%d_%s_%s_%s_%s_s.png"%(w0,m0,r,s,mu,t),dpi=300,bbox_inches='tight')
	plt.show()
	'''

	#print(s_max)
	from extreme_value_statistics import *
	
	#fitting to EVT	
	s_ran=np.arange(np.min(s_max),np.max(s_max),0.1)
	
	

	ax=plt.axes((0.1,0.1,0.4,0.3))
	smmean=np.mean(s_max)
	smstd=np.std(s_max)
	#print(smmean,smstd)
	s_st=(s_max-smmean)/smstd
	ax.hist(s_st,density=True,bins=np.arange(np.min(s_st),np.max(s_st),0.01)*1.1,alpha=0.8)
	x=np.arange(-2,10,0.01)
	ax.plot(x,Gumbel_pdf(x),label='Gumbel',c='red')
	ax.plot(x,Frechet_pdf(x,gamma=1.01),label='Frechet',c='black')
	ax.set_xlabel(r'$s_\mathrm{max}$')
	ax.set_ylabel(r'$P(s_\mathrm{max})$')
	ax.legend(frameon=False)
	#ax.set_yscale('log')
	#ax.set_ylim(ymin=10**-7)
	#plt.savefig("%d_%d_%s_%s_%s_%s_s_max.png"%(w0,m0,r,s,mu,t),dpi=300,bbox_inches='tight')
	plt.show()	
	'''
