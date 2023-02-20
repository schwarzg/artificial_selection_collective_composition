import numpy as np
import matplotlib.pyplot as plt
import mpltern

#parameter prepare
mu=1e-4
r=0.5
s=3e-2
N0=1000
mbar=100
ncomm=10
g=ncomm
rhat=0
nensemble=40
ncycle=10
tcycle=np.log(ncomm+1)/r

#######################
# maturation flow and accessible region
#######################

ax=plt.axes((0.8,0.7,0.30,0.35),projection='ternary')

ax.annotate('b',xy=(-0.3,1.35),weight='bold',xycoords='axes fraction')
ax.annotate(' Composition flow in maturation step and accessible region',xy=(-0.25,1.35),xycoords='axes fraction')
def dfdt_vec(p):
	#print('p',p)
	z,x,y=p #fv,fw,fm
	#x=x/scale
	#y=y/scale
	#z=z/scale
	m=s*y*(1-y)-2*s*y*z+mu*(1-2*y-z)
	v=2*s*z*(1-z)-s*y*z+mu*y
	w=-m-v

	return np.array([v,w,m])


from mpltern.ternary.datasets import get_triangular_grid,get_shanon_entropies
_,_,_,v=get_shanon_entropies()
t,l,R=get_triangular_grid()
dt,dl,dr=dfdt_vec((t,l,R))

length=np.sqrt(dl**2+dr**2+dt**2)

data=np.loadtxt("N0%d_r%s_s%s_mu%s_ext_means.dat"%(N0,r,s,mu))
scale=1
fms=data[:,0]
fvs=data[:,1]
fws=1-fms-fvs
v1x=data[:,2]-data[:,0]
v1y=data[:,3]-data[:,1]
v2x=data[:,4]-data[:,0]
v2y=data[:,5]-data[:,1]

mask=fws>-1e-8
fws=fws[mask]
fms=fms[mask]
fvs=fvs[mask]
v1x=v1x[mask]
v2x=v2x[mask]
v1y=v1y[mask]
v2y=v2y[mask]
dscore=np.zeros(len(v1x))
for i,fm in enumerate(fms):
	midx=int(fm)
	vidx=int(fvs[i])
	#dscore[(midx,vidx)]=score[i]
	if np.sign(v1x[i])*np.sign(v2x[i]-np.sign(v1x[i])*mu*fms[i])<=0 and np.sign(v1y[i])*np.sign(v2y[i]-np.sign(v1y[i])*mu*fvs[i])<=0:
	#if np.sign(v1x[i])*np.sign(v2x[i])<=0 and np.sign(v1y[i])*np.sign(v2y[i])<=0:
	#if np.sign(v1y[i])*np.sign(v2y[i])<0:
		dscore[i]=1
	else:
		dscore[i]=0
import matplotlib as mpl
oldgrey=mpl.cm.get_cmap('Greys')
newgrey=mpl.colors.ListedColormap(oldgrey(np.linspace(0,0.5,3)))
tc=ax.tripcolor(fvs,fws,fms,dscore,shading='gouraud',rasterized=True,cmap=newgrey)

qv=ax.quiver(t,l,R,dt,dl,dr,length)
cax2 = ax.inset_axes([1.3, 0.1, 0.05, 0.9], transform=ax.transAxes)
cbar2=plt.colorbar(qv,cax=cax2)
cbar2.set_label('Speed',rotation=270,va='baseline')

ax.set_tlabel('Double-mutant frequency')
ax.set_llabel('Wild-type frequency')
ax.set_rlabel('Single-mutant frequency')

ax.taxis.set_ticks_position('tick1')
ax.laxis.set_ticks_position('tick1')
ax.raxis.set_ticks_position('tick1')

ax.taxis.set_label_position('tick1')
ax.laxis.set_label_position('tick1')
ax.raxis.set_label_position('tick1')
#######################
# Schematic explanation of how to get the region
#######################

################# preparation
cc=['#75140c','#0044aa','#a004bf']

x1=[0,1]
y1=[0,0]
x2=[0,0.5]
y2=[0,np.sqrt(3)*0.5]
x3=[0.5,1]
y3=[np.sqrt(3)*0.5,0]

cv1yl=0.05
cv1xl=cv1yl/np.sqrt(3)
cv1yr=0.25
cv1xr=1-cv1yr/np.sqrt(3)

def curve1(x,a,b):
	return a*x**2+b
def curve1p(x,x1,y1,x2,y2):
	return (y2-y1)/(x2-x1)**4*(x-x1)**4+y1

a1=(cv1yl-cv1yr)/(cv1xl**2-cv1xr**2)
b1=(-cv1xr**2*cv1yl+cv1xl**2*cv1yr)/(cv1xl**2-cv1xr**2)

cv1x=np.arange(cv1xl,cv1xr+0.005,0.01)
cv1y=curve1p(cv1x,cv1xl,cv1yl,cv1xr,cv1yr)

cv2yl=0.61
cv2xl=cv2yl/np.sqrt(3)
cv2yr=0.6
cv2xr=1-cv2yr/np.sqrt(3)

def curve2(x,a,b,c):
	return a*(x-c)**2+b

c=0.51
a2=(cv2yl-cv2yr)/((cv2xl-c)**2-(cv2xr-c)**2)
b2=(-(cv2xr-c)**2*cv2yl+(cv2xl-c)**2*cv2yr)/((cv2xl-c)**2-(cv2xr-c)**2)

cv2x=np.arange(cv2xl,cv2xr+0.005,0.01)
cv2y=curve2(cv2x,a2,b2,c)

#low region1
xl1=np.linspace(0,cv1xl)
yl1=np.sqrt(3)*xl1

#low region2
xl2=np.linspace(cv1xl,cv1xr)
#yl2=curve1(xl2,a1,b1)
yl2=curve1p(xl2,cv1xl,cv1yl,cv1xr,cv1yr)

#low region3
xl3=np.linspace(cv1xr,1)
yl3=-np.sqrt(3)*(xl3-1)


#high region1
xh1=np.linspace(cv2xl,0.5)
yh1u=np.sqrt(3)*xh1
yh1d=curve2(xh1,a2,b2,c)


#high region1
xh2=np.linspace(0.5,cv2xr)
yh2u=-np.sqrt(3)*(xh2-1)
yh2d=curve2(xh2,a2,b2,c)


#########################basic location
bxx=0.05
bxy=0.70

#square frame
bxsq=plt.axes((bxx-0.03,bxy-0.28,0.63,0.75))
bxsq.tick_params(left=False,bottom=False,labelleft=False,labelbottom=False)
bxsq.annotate('a',xy=(-0.0,1.05),weight='bold',xycoords='axes fraction')
bxsq.annotate(' Prediction of accessible region by collective selection',xy=(0.05,1.05),xycoords='axes fraction')

def draw_filled_circle(ax,cenx,ceny,radx,rady=None,center='on',angle=0,label=''):
	#draw unit circle
	
	t=np.linspace(-np.pi,0,20)
	if rady is None:
		rady=radx
	
	xu=radx*np.cos(t)
	xd=radx*np.cos(t)
	yu=-rady*np.sin(t)
	yd=rady*np.sin(t)

	if angle!=0:
		#rotation, unit : pi
		#angle=angle%np.pi
		s=np.sin(angle)
		c=np.cos(angle)
		rot_xu=c*xu-s*yu
		rot_yu=s*xu+c*yu
		rot_xd=c*xd-s*yd
		rot_yd=s*xd+c*yd
		xu=rot_xu	
		xd=rot_xd	
		yu=rot_yu	
		yd=rot_yd	

	
	#translation
	yubot=xu*np.tan(angle)+ceny
	ydtop=xd*np.tan(angle)+ceny
	xu=xu+cenx
	xd=xd+cenx
	yu=yu+ceny
	yd=yd+ceny


	pc=ax.fill_between(xu,yu,yubot,alpha=0.5)
	ax.fill_between(xd,ydtop,yd,fc=pc.get_fc(),alpha=0.5,label=label)
	if center=='on':
		ax.scatter(cenx,ceny,c=pc.get_fc(),alpha=1)
	return ax

bx0=plt.axes((bxx-0.05,bxy+0.06,0.35,0.35))
bx1=plt.axes((bxx+0.28,bxy+0.06,0.35,0.35))
bx2=plt.axes((bxx,bxy-0.28,0.25,0.25))
bx3=plt.axes((bxx+0.32,bxy-0.28,0.25,0.25))

bx0.plot(x1,y1,c='black')
bx0.plot(x2,y2,c='black')
bx0.plot(x3,y3,c='black')
bx0.add_patch(mpl.patches.Circle((-0.10,-0.10), radius=0.07, color=cc[0],ec='black'))
bx0.add_patch(mpl.patches.Circle((1.10,-0.10), radius=0.07, color=cc[1],ec='black'))
bx0.add_patch(mpl.patches.Circle((0.5,1.00), radius=0.07, color=cc[2],ec='black'))
bx0.scatter(0.2,0.1,marker='o',fc=None)
bx0.axis('scaled')
bx0.annotate('',xy=(0.45,0.25),xytext=(0.2,0.1),arrowprops=dict(arrowstyle='->'))
bx0=draw_filled_circle(bx0,0.45,0.25,0.25,0.15,angle=np.pi/6.0)
bx0.annotate('1. Adult distribution',xy=(0,1.05),xycoords='axes fraction')
bx0.text(
    1.1, 0.5, "Select", ha="left", va="center", size=10,
    bbox=dict(boxstyle="rarrow,pad=0.3", fc='w',lw=1))
bx0.axis('off')

bx1.plot(x1,y1,c='black')
bx1.plot(x2,y2,c='black')
bx1.plot(x3,y3,c='black')
bx1.add_patch(mpl.patches.Circle((-0.10,-0.10), radius=0.07, color=cc[0],ec='black'))
bx1.add_patch(mpl.patches.Circle((1.10,-0.10), radius=0.07, color=cc[1],ec='black'))
bx1.add_patch(mpl.patches.Circle((0.5,1.00), radius=0.07, color=cc[2],ec='black'))

bx1.scatter(0.2,0.1,marker='o',fc=None,label='Initial')
bx1=draw_filled_circle(bx1,0.45,0.25,0.25,0.15,angle=np.pi/6.0,label='Adult')
bx1=draw_filled_circle(bx1,0.30,0.16,0.08,0.10,angle=np.pi/6.,label='Selected')

bx1.axis('scaled')
bx1.annotate('2. Selected distribution',xy=(0,1.05),xycoords='axes fraction')
bx1.legend(loc=(-1,-0.1),ncol=3,frameon=False)
bx1.axis('off')

bx2.scatter(0.36,0.2,marker='o',zorder=30,fc=None)
bx2=draw_filled_circle(bx2,0.45,0.25,0.25,0.15,angle=np.pi/6.0)
bx2=draw_filled_circle(bx2,0.30,0.16,0.08,0.10,angle=np.pi/6.)
bx2.annotate('',xy=(0.45,0.25),xytext=(0.36,0.2),arrowprops=dict(arrowstyle='->'))
bx2.annotate('',xy=(0.30,0.16),xytext=(0.36,0.2),arrowprops=dict(arrowstyle='->'))
bx2.axis('scaled')
bx2.annotate('3. Case: Accessable ',xy=(-0.1,1.05),xycoords='axes fraction')
bx2.axis('off')

bx3.scatter(0.25,0.1,marker='o',zorder=30,fc=None)
bx3=draw_filled_circle(bx3,0.45,0.25,0.25,0.15,angle=np.pi/6.0)
bx3=draw_filled_circle(bx3,0.30,0.16,0.08,0.10,angle=np.pi/6.)
bx3.annotate('',xy=(0.45,0.25),xytext=(0.25,0.1),arrowprops=dict(arrowstyle='->'))
bx3.annotate('',xy=(0.30,0.16),xytext=(0.25,0.1),arrowprops=dict(arrowstyle='->'))
bx3.axis('scaled')
bx3.annotate('4. Case: Not accessible ',xy=(-0.1,1.05),xycoords='axes fraction')
bx3.axis('off')
'''
bx3.annotate('4. Compare directions',xy=(0,1.05),xycoords='axes fraction')
bx3.annotate('',xy=(3,2),xytext=(2,2),arrowprops=dict(arrowstyle='->'))
bx3.annotate('',xy=(1,2.1),xytext=(2,2),arrowprops=dict(arrowstyle='->'))
bx3.annotate('',xy=(2,0.7),xytext=(1,0.5),arrowprops=dict(arrowstyle='->'))
bx3.annotate('',xy=(3,0.3),xytext=(1,0.5),arrowprops=dict(arrowstyle='->'))
bx3.set_xlim(0,4)
bx3.set_ylim(0,3)
bx3.scatter([2,1],[2,0.5]) #initial
bx3.scatter([3,3],[2,0.3]) #maturation
bx3.scatter([1,2],[2.1,0.7]) #selection
bx3.annotate('Opposite : Accessible',xy=(2,2.5),ha='center',va='bottom')
bx3.annotate('In phase : Not accessible',xy=(2,1),ha='center',va='bottom')
bx3.axis('off')
'''

#############################
# Plots of numerical tests
#############################
import itertools as itt
def cf0(w,m,v):
	#community function for m 
	N=w+m+v
	f=np.divide(w,N)
	return f
def cf1(w,m,v):
	#community function for m 
	N=w+m+v
	f=np.divide(m,N)
	return f
def cf2(w,m,v):
	#community function for v
	N=w+m+v
	f=np.divide(v,N)
	return f
def Draw_triangle_from_data(ax,m0,v0,mhats,vhats):
	#draw background
	ax.tripcolor(fvs,fws,fms,dscore,shading='gouraud',rasterized=True,cmap=newgrey)

	for (mhat,vhat) in itt.product(mhats,vhats):
	
		#Initia/Target point
		fm=m0/N0
		fv=v0/N0
		fw=1-fm-fv
		ax.scatter(fv,fw,fm,marker='o',s=20,ec='black',fc='white',zorder=10,label='Initial')
		ax.scatter(vhat,1-vhat-mhat,mhat,marker='o',s=20,c='black',zorder=10,label='Target')
	
		AGS=[]
		folder="data/raw/"
		nensemble=30
		for e in range(nensemble):
			AGS.append(np.loadtxt(folder+"N0%s_m0%s_v0%s_r%s_s%s_mu%s_g%s_mhat%s_vhat%s_AGS_point_%d.cycle"%(N0,m0,v0,r,s,mu,ncomm,mhat,vhat,e)))
		AGS=np.array(AGS) # [ensemble, timestamp, (T,selected index,w----,m----,v----)]
		
		T=AGS[0,:,0]
		sel_inds=AGS[:,:,1].astype(int) #selected index - [ensemble, timestamp]
		c1_sel=np.zeros((nensemble,len(T),ncomm)) #c avg datas, [ensemble, timestamp]
		c2_sel=np.zeros((nensemble,len(T),ncomm)) #c avg datas, [ensemble, timestamp]
		w_cho=np.zeros((nensemble,len(T))) #w datas, [ensemble, timestamp]
		m_cho=np.zeros((nensemble,len(T))) #m datas, [ensemble, timestamp]
		v_cho=np.zeros((nensemble,len(T))) #m datas, [ensemble, timestamp]
		#ensemble
		for i in range(nensemble): #timestamp
			for j in range(len(T)):
				#c_sel[i,j,:]=cf(AGS[i,j,ncomm+2:ncomm*2+2],AGS[i,j,ncomm*3+2:])
				w_cho[i,j]=AGS[i,j,sel_inds[i,j]+2+ncomm]
				m_cho[i,j]=AGS[i,j,sel_inds[i,j]+2+ncomm*3]
				v_cho[i,j]=AGS[i,j,sel_inds[i,j]+2+ncomm*5]
		
		#Artificial group selection with the selected
		c0_cho=cf0(w_cho,m_cho,v_cho)
		c1_cho=cf1(w_cho,m_cho,v_cho)
		c2_cho=cf2(w_cho,m_cho,v_cho)
		cho0avg=np.mean(c0_cho,axis=0)	
		cho1avg=np.mean(c1_cho,axis=0)	
		cho2avg=np.mean(c2_cho,axis=0)
		ax.plot(cho2avg,cho0avg,cho1avg)	

	return ax

#cxsq=plt.axes((0.1,-1.0,1.4,1.6))

cx11=plt.axes((0.16,0.0,0.20,0.23),projection='ternary')
cx11.annotate('c',xy=(-0.1,1.05),xycoords='axes fraction',weight='bold')
cx11=Draw_triangle_from_data(cx11,85,5,[0.7],[0.15])
cx11.legend(frameon=False,loc=(1.0,1.0))

cx12=plt.axes((0.56,0.0,0.20,0.23),projection='ternary')
cx12.annotate('d',xy=(-0.1,1.05),xycoords='axes fraction',weight='bold')
cx12=Draw_triangle_from_data(cx12,85,5,[0.33],[0.33])

cx13=plt.axes((0.94,0.0,0.20,0.23),projection='ternary')
cx13.annotate('e',xy=(-0.1,1.05),xycoords='axes fraction',weight='bold')
cx13=Draw_triangle_from_data(cx13,85,5,[0.02],[0.9])

cx21=plt.axes((0.16,-0.40,0.20,0.23),projection='ternary')
cx21.annotate('f',xy=(-0.1,1.05),xycoords='axes fraction',weight='bold')
cx21=Draw_triangle_from_data(cx21,50,50,[0.7],[0.15])

cx22=plt.axes((0.56,-0.40,0.20,0.23),projection='ternary')
cx22.annotate('g',xy=(-0.1,1.05),xycoords='axes fraction',weight='bold')
cx22=Draw_triangle_from_data(cx22,50,50,[0.33],[0.33])

cx23=plt.axes((0.94,-0.40,0.20,0.23),projection='ternary')
cx23.annotate('h',xy=(-0.1,1.05),xycoords='axes fraction',weight='bold')
cx23=Draw_triangle_from_data(cx23,50,50,[0.02],[0.9])

cx31=plt.axes((0.16,-0.80,0.20,0.23),projection='ternary')
cx31.annotate('i',xy=(-0.1,1.05),xycoords='axes fraction',weight='bold')
cx31=Draw_triangle_from_data(cx31,50,900,[0.7],[0.15])

cx32=plt.axes((0.56,-0.80,0.20,0.23),projection='ternary')
cx32.annotate('j',xy=(-0.1,1.05),xycoords='axes fraction',weight='bold')
cx32=Draw_triangle_from_data(cx32,50,900,[0.33],[0.33])

cx33=plt.axes((0.94,-0.80,0.20,0.23),projection='ternary')
cx33.annotate('k',xy=(-0.1,1.05),xycoords='axes fraction',weight='bold')
cx33=Draw_triangle_from_data(cx33,50,900,[0.02],[0.9])

formatter='svg'
plt.savefig('figures/FigS6.'+formatter,bbox_inches='tight',dpi=300,format=formatter)
plt.show()

