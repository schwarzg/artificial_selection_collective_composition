import numpy as np

def barw_th(w0,t,r,mu):
	return w0*np.exp((r-mu)*t)	
def barm_th(m0,w0,t,r,mu,s):
	erst=np.exp((r+s)*t)
	return m0*erst+mu*w0*(erst-np.exp((r-mu)*t))/(s+mu)

#Insert variance -> return variance
def sig2w_th(w0,sig2w0,t,r,mu):
	ermt=np.exp((r-mu)*t)
	return sig2w0*ermt**2+(r+mu)*w0*ermt*(ermt-1)/(r-mu)

#Insert variance -> return variance insert sigwm0, basically same with -sig2m0 if w0+m0=N0
def sig2m_th(m0,sig2m0,w0,sig2w0,covwm0,t,r,mu,s):
	erst=np.exp((r+s)*t)
	erst2=erst**2
	ermt=np.exp((r-mu)*t)
	O1=sig2m0*erst2
	Omu=(m0+mu*w0/(s+mu))*erst*(erst-1)-(r-mu)*mu*w0*(erst2-ermt)/(r+2*s+mu)/(s+mu)+2*mu*covwm0*erst*(erst-ermt)/(s+mu)
	Omu2=mu**2*(sig2w0+(r+mu)*w0/(r-mu))*(erst-ermt)**2/(s+mu)**2
	Omu3=-4*mu**2*r*w0*((erst2-erst)/(s+mu)-(erst-ermt)/(r+2*s+mu))/(r-mu)/(r+s)
	return O1+Omu+Omu2+Omu3

#Gaussian approximations
def barc_th(w0,m0,t,r,mu,s):
	w=barw_th(w0,t,r,mu)
	m=barm_th(m0,w0,t,r,mu,s)
	N=w+m
	return m/N
def barc_th_v2(f0,N0,t,r,mu,s):
	w=barw_th((1-f0)*N0,t,r,mu)
	m=barm_th(f0*N0,(1-f0)*N0,t,r,mu,s)
	N=w+m
	return m/N
def barc_th_v3(f0,t,r,mu,s):
	w=(1-f0)
	est=np.exp(s*t)
	m=f0*est+mu/s*(1-f0)*(est-1)
	return m/(w+m)

#initial 'variances' -> return variance
def sig2c_th(w0,sig2w0,m0,sig2m0,covwm0,t,r,mu,s):
	w=barw_th(w0,t,r,mu)
	m=barm_th(m0,w0,t,r,s,mu)
	N=(w+m)
	N2=N**2
	c=m/N
	sig2w=sig2w_th(w0,sig2w0,t,r,mu)
	sig2m=sig2m_th(m0,sig2m0,w0,sig2w0,covwm0,t,r,mu,s)
	return ((1-c)**2*sig2m+c**2*sig2w)/N2	

#Assume sig2w=sig2m at intial 'variances' -> return variance
def sig2c_th_v2(f0,sig2f0,N0,t,r,mu,s):

	#print(f0,sig2f0,N0,t,r,mu,s)
	m0=f0*N0
	w0=(1-f0)*N0
	sig2w0=(N0**2)*sig2f0/(2*f0**2-2*f0+1)
	sig2m0=sig2w0
	#print(w0,m0,sig2w0,sig2m0)
	
	w=barw_th(w0,t,r,mu)
	m=barm_th(m0,w0,t,r,s,mu)
	N=(w+m)
	N2=N**2
	c=m/N
	sig2w=sig2w_th(w0,sig2w0,t,r,mu)
	sig2m=sig2m_th(m0,sig2m0,w0,sig2m0,-sig2m0,t,r,mu,s)
	#print(w,m,sig2w,sig2m)
	return ((1-c)**2*sig2m+c**2*sig2w)/N2	

#nascent dirac delta function in probabilistic sence
def ndirac(x,eps=1e-4):
		return 1./eps if np.fabs(x)<eps/2 else 0

#Q(1/g) from given gaussian distribution (f0,sig0)
def freq_min_th(g,f0,sig0):
	return st.norm.ppf(1/g,loc=f0,scale=sig0)	

#Reproduction variance
def barRm_th(c,N0=1000):
	return np.array(N0*c).astype(int)	

#return population 'variance'	
def sigRm_th(c,barN,N0=1000):
	return N0*c*(1.-c)*(barN-N0)/(barN-1)

#return sample 'variance'
def sigRmbar_th(c,barN,N0,g):
	return N0*c*(1.-c)*(barN-g*N0)/(barN-1)/g


if __name__=='__main__':
	import matplotlib.pyplot as plt
	C=np.linspace(0,1,100)
	barN=11*1000
	plt.plot(C,sigRm_th(C,barN))
	plt.show()

	T=np.arange(0,5,0.1)
	w0=800
	m0=200
	sig2w0=m0
	sig2m0=m0
	f=m0/(m0+w0)
	N=(m0+w0)
	sig0=((1-f)**2*sig2m0+f**2*sig2w0)/N**2	
	print('sigf',sig0)
	sigcs=sig2c_th(w0,sig2w0,m0,sig2m0,-sig2m0,T,0.5,0.0001,0.03)
	sigcs2=sig2c_th_v2(f,sig0,N,T,0.5,0.0001,0.03)
	plt.plot(T,sigcs)
	plt.plot(T,sigcs2)
	plt.show()

	from model import model
	from tau_leaping import tau_leaping
	from growth_alter import *
	import scipy.stats as st
	
	mu=1e-4
	r=0.5
	s=-3e-2
	N0=1000
	mbar=100
	ncomm=10
	rhat=0.05
	nens=300
	ncycle=10
	tcycle=np.log(ncomm+1)/r

	proFunc=[
	lambda x: r*x[0],
	lambda x: mu*x[0],
	lambda x: (r+s)*x[1]
	]
	changeVec=np.array([[1,0],[-1,1],[0,1]]).astype(int)
	rxnOrder=np.array([[1,0],[1,0],[0,1]]).astype(int)

	growthmodel=model(proFunc=proFunc,changeVec=changeVec,rxnOrder=rxnOrder)
	
	solver=tau_leaping(model=growthmodel,eps=0.01)

	w0=800
	m0=200

	m0=np.ones(nens)*800#np.clip(np.random.poisson(200,size=nens),0,1000)
	w0=1000-m0
	
	barw=barw_th(np.mean(w0),tcycle,r,mu)	
	sig2w=sig2w_th(np.mean(w0),np.var(w0),tcycle,r,mu)
	barm=barm_th(np.mean(m0),np.mean(w0),tcycle,r,mu,s)
	sig2m=sig2m_th(np.mean(m0),np.var(m0),np.mean(w0),np.var(w0),-np.var(m0),tcycle,r,mu,s)	
	def get_f(w,m):
		return np.divide(m,w+m)
	c0=get_f(w0,m0)
	barc=barc_th_v2(np.mean(c0),N0,tcycle,r,mu,s)
	sig2c=sig2c_th_v2(np.mean(c0),np.var(c0),N0,tcycle,r,mu,s)
	
	ws=np.zeros(nens)
	ms=np.zeros(nens)
	ws2=np.zeros(nens)
	ms2=np.zeros(nens)
	for e in range(nens):
		print(e)
		_,X=solver.run(np.array([w0[e],m0[e]]),0,tcycle)
		ws[e]=X[0,-1]	
		ms[e]=X[1,-1]	
		ws2[e]=growth_sampling_w(np.mean(w0),np.var(w0),tcycle,r,mu,s)
		ms2[e]=growth_sampling_m(np.mean(m0),np.var(m0),np.mean(w0),np.var(w0),-np.var(m0),tcycle,r,mu,s)
	cs=get_f(ws,ms)	
	#ax=plt.axes((0.1,0.5,0.4,0.3))
	ax2=plt.axes((0.1,0.5,0.35,0.3))
	#ax.hist(w0,bins=10,density=True)
	ax2.hist(ws,bins=30,density=True,alpha=0.4,label='tau')
	ax2.hist(ws2,bins=30,density=True,alpha=0.4,label='samp')
	wrange=np.linspace(barw-3*np.sqrt(sig2w),barw+3*np.sqrt(sig2w),100)	
	print(barw,sig2w,'?=',np.mean(ws2),np.var(ws2))	
	pdfw=st.norm.pdf(wrange,loc=barw,scale=np.sqrt(sig2w))
	ax2.plot(wrange,pdfw,label='Gauss')
	ax2.set_xlabel(r'$w$')
	ax2.set_ylabel(r'$P(w)$')
	ax2.legend(frameon=False,handlelength=0.5)

	#bx=plt.axes((0.1,0.1,0.4,0.3))
	bx2=plt.axes((0.5,0.5,0.35,0.3))
	#bx.hist(m0,bins=10,density=True)
	bx2.hist(ms,bins=30,density=True,alpha=0.4)		
	bx2.hist(ms2,bins=30,density=True,alpha=0.4)		
	mrange=np.linspace(barm-3*np.sqrt(sig2m),barm+3*np.sqrt(sig2m),100)	
	print(barm,sig2m,'?=',np.mean(ms2),np.var(ms2))	
	pdfm=st.norm.pdf(mrange,loc=barm,scale=np.sqrt(sig2m))
	bx2.plot(mrange,pdfm)
	bx2.set_xlabel(r'$m$')
	bx2.set_ylabel(r'$P(m)$')
	cx2=plt.axes((0.1,0.1,0.35,0.3))
	cx2.hist(cs,density=True,alpha=0.4)
	crange=np.linspace(barc-3*np.sqrt(sig2c),barc+3*np.sqrt(sig2c),100)	
	print(barc,sig2c,'?=',np.mean(cs),np.var(cs))	
	pdfc=st.norm.pdf(crange,loc=barc,scale=np.sqrt(sig2c))
	cx2.plot(crange,pdfc)
	#plt.savefig('compare_tausampgauss.svg',dpi=300,bbox_inches='tight',format='svg')
	plt.show()
