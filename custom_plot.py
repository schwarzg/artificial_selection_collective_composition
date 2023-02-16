import matplotlib.pyplot as plo
import numpy as np

def violinplot_from_histogram(ax,hists,bins,positions=None,showmeans=True,width=0.5,side='right',yoffset=None,marker='o',**kwargs):

	nhists=len(hists)
	if type(positions) is type(None):
		positions=range(nhists) 
	if type(yoffset) is type(None):
		yoffset=np.zeros(nhists)

	#General : bins for each histograms are given	
	if len(bins.shape)!=1:
		dx=bins[:,1]-bins[:,0]
		x=bins[:,:-1]+0.5*np.array([dx]).T
	else:
		dx=bins[1]-bins[0]
		x=bins[:-1]+0.5*dx
		x=np.repeat([x],nhists,axis=0)
		dx=np.ones(nhists)*dx
	means=[]
	for ind,dats in enumerate(hists):
		#Normalize for safe
		Norm=np.sum(dats*dx[ind])
		dats=dats/Norm
		maxval=np.max(dats)
		#check meaningful range
		#effrange=range(len(x))
		if side=='right' or side=='both':
			ax.fill_between(positions[ind]+0.5*dats/maxval*width,0,x[ind]-yoffset[ind],**kwargs)
		if side=='left' or side=='both':
			ax.fill_between(positions[ind]-0.5*dats/maxval*width,0,x[ind]-yoffset[ind],**kwargs)
		if showmeans==True:
			meanval=np.sum(x[ind]*dats*dx[ind])
			print(meanval)
			means.append(meanval)
			ax.scatter(positions[ind],meanval-yoffset[ind],**kwargs,marker=marker)
	return ax,means	

import scipy.stats as st
def statistical_significant(ax,dat1,dat2,pos,offset=0,upside=True,**kwargs):
	#pos=[dat1,dat2]
	
	#calulate p value
	_,pv=st.ttest_ind(dat1,dat2,equal_var=False)
	if pv<0.001:
		nstar=3
	elif pv<0.01:
		nsatr=2	
	elif pv<0.05:
		nstar=1
	else:
		nstar=0

	if upside==True:	
		#Vertical lines positions
		m1=np.max(dat1)
		m2=np.max(dat2)
		by=np.maximum(m1,m2)
		disp=ax.transAxes.transform((1,1))
		print(disp)
		dy=ax.transData.inverted().transform(disp)
		print(dy)
		ty=by+np.abs(0.1*by)#+dy[1]
		
		#print(m1,m2,by,ty)
	
		#calculate height of horizontal bar
		annotx=0.5*(pos[0]+pos[1])
		annoty=ty
	
		#draw frames
		ax.vlines(pos,by,ty,colors='black')
		ax.hlines(ty,pos[0],pos[1],colors='black')
		
		txt='*'*nstar
		if nstar==0:
			txt='NS'
		print(txt)
		ax.annotate(txt,(annotx,annoty),horizontalalignment='center')	
	
	else: # Downside		
		#Vertical lines positions
		m1=np.min(dat1)
		m2=np.min(dat2)
		ty=np.minimum(m1,m2)
		disp=ax.transAxes.transform((1,1))
		print(disp)
		dy=ax.transData.inverted().transform(disp)
		print(dy)
		by=ty-np.abs(0.1*ty)#+dy[1]
		
		#print(m1,m2,by,ty)
	
		#calculate height of horizontal bar
		annotx=0.5*(pos[0]+pos[1])
		annoty=by
	
		#draw frames
		ax.vlines(pos,by,ty,colors='black')
		ax.hlines(by,pos[0],pos[1],colors='black')
		
		txt='*'*nstar
		if nstar==0:
			txt='NS'
		print(txt)
		ax.annotate(txt,(annotx,annoty),horizontalalignment='center',verticalalignment='bottom')
		

	return ax
