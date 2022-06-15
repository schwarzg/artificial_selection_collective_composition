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
