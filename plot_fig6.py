import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

############################################
# Schematic - maturation of three strain
############################################

cc=['#75140c','#0044aa','#a004bf']

#location fix
mxx=0.1
mxy=0.45-0.28-0.05
mx=plt.axes((mxx,mxy,0.45,0.45))
mx.annotate('a',xy=(0.,1.075),weight='bold',xycoords='axes fraction')
mx.annotate(' Three-strain system',xy=(0.05,1.075),xycoords='axes fraction')
def draw_birth(ax,x,y,ut='',dt='',c='r',scale=1,length=3):
    rad=0.4
    ax.add_patch(mpl.patches.Circle((x,y), radius=rad, color=c,ec='black'))
    ax.annotate('',xy=(x+0.5+length,y),xytext=(x+0.5,y),arrowprops=dict(arrowstyle='->'))
    ax.add_patch(mpl.patches.Circle((x+length+1,y), radius=rad, color=c,ec='black'))
    ax.add_patch(mpl.patches.Circle((x+length+2,y), radius=rad, color=c,ec='black'))
    ax.annotate(ut,xy=(x+2,y+0.1),ha='center',va='bottom')
    ax.annotate(dt,xy=(x+2,y-0.1),ha='center',va='top')
    return ax

def draw_mutation(ax,x,y,cl='r',cr='b',scale=1):
    rad=0.4
    ax.add_patch(mpl.patches.Circle((x,y), radius=rad, color=cl,ec='black'))
    ax.annotate('',xy=(x+3.5,y),xytext=(x+0.5,y),arrowprops=dict(arrowstyle='->'))
    ax.add_patch(mpl.patches.Circle((x+4,y), radius=rad, color=cr,ec='black'))
    ax.annotate('',xy=(x+2,y+0.1),ha='center',va='bottom')
    return ax

rx1=0
ry1=0

mx.annotate(r'Reactions',xy=(rx1,ry1+2.),ha='left',va='center')
mx.annotate(r'Rates',xy=(rx1+6.5,ry1+2.),ha='left',va='center')


draw_birth(mx,rx1,ry1,c=cc[0],length=3)
mx.annotate(r'$r$',xy=(rx1+6.5,ry1),ha='left',va='center')


ry2=-2.5
draw_birth(mx,rx1,ry2,c=cc[1],length=3)
mx.annotate(r'$r+s$',xy=(rx1+6.5,ry2),ha='left',va='center')

ry3=-5
draw_birth(mx,rx1,ry3,c=cc[2],length=3)
mx.annotate(r'$r+2s$',xy=(rx1+6.5,ry3),ha='left',va='center')

rx2=7
ry4=ry3-3.0
draw_mutation(mx,rx1,ry4,cl=cc[0],cr=cc[1])
mx.annotate(r'$\mu$',xy=(rx1+6.5,ry4),ha='left',va='center')

ry5=ry3-5.5
draw_mutation(mx,rx1,ry5,cl=cc[1],cr=cc[2])
mx.annotate(r'$\mu$',xy=(rx1+6.5,ry5),ha='left',va='center')

mx.axis('scaled')
mx.tick_params(left=False,bottom=False,labelleft=False,labelbottom=False)
mx.set_ylim(ymin=-12.5,ymax=3)
mx.set_xlim(xmin=-2,xmax=10.0)

    

    

############################################
#prepare for triangle plot
############################################

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

#point initial and target
i1x=0.5
i1y=0.75
i2x=0.05
i2y=0.03

i3x=0.3
i3y=0.00


#Three strain case
#cxx=0.1
#cxy=-0.30
cxx=mxx+0.35#0.55
cxy=mxy+0.1#0.37

cx=plt.axes((cxx,cxy,0.35,0.35))
cx.annotate('b',xy=(-0.0,1.1),weight='bold',xycoords='axes fraction')
cx.annotate(' Simulation of composition trajectories according to initial/target points',xy=(0.05,1.1),xycoords='axes fraction')
cx.plot(x1,y1,c='black')
cx.plot(x2,y2,c='black')
cx.plot(x3,y3,c='black')
cx.add_patch(mpl.patches.Circle((-0.10,-0.10), radius=0.07, color=cc[0],ec='black'))
cx.add_patch(mpl.patches.Circle((1.10,-0.10), radius=0.07, color=cc[1],ec='black'))
cx.add_patch(mpl.patches.Circle((0.5,1.00), radius=0.07, color=cc[2],ec='black'))


cx.fill_between(xl1,yl1,0,color='lightgray')
cx.fill_between(xl2,yl2,0,color='lightgray')
cx.fill_between(xl3,yl3,0,color='lightgray')
cx.fill_between(xh1,yh1u,yh1d,color='lightgray')
cx.fill_between(xh2,yh2u,yh2d,color='lightgray',label='Accessible by inter-collective selection')
cx.fill_between(xh2,yh2u,yh2u,color='white',label='Not accessible')

cx.scatter(i1x,i1y,color='black',label='Initial',fc='white',zorder=20)
cx.scatter(i2x,i2y,color='black',fc='white',zorder=20)
cx.scatter(i3x,i3y,color='black',fc='white',zorder=20)
fx=0.4
fy=np.sqrt(3)*0.37
cx.scatter(fx,fy,color='black',label='Target')

cx.annotate('',xy=(fx,fy),xytext=(i1x,i1y),arrowprops=dict(arrowstyle='->',color='green'),label='')
cx.annotate('',xy=(fx,fy),xytext=(i2x,i2y),arrowprops=dict(arrowstyle='->',color='green'))
cx.annotate('',xy=(fx,fy),xytext=(i3x,i3y),arrowprops=dict(arrowstyle='->',color='green'))
lb=cx.legend(frameon=False,loc=(0,-0.3),ncol=2)
lb.legendHandles[1].set_edgecolor('black')
cx.axis('scaled')
cx.axis('off')

cx.annotate('',xy=(0.8,0.1),xytext=(0.2,0.1),xycoords='axes fraction',arrowprops=dict(arrowstyle='simple',color='C0',mutation_scale=20))
cx.annotate('',xy=(0.9,0.8),xytext=(0.9,0.2),xycoords='axes fraction',arrowprops=dict(arrowstyle='simple',color='C0',mutation_scale=20))
cx.annotate('Intra-collective\nSelection',xy=(1.1,0.9),ha='center',va='bottom',color='C0')


#second plot
cx2=plt.axes((cxx+0.30,cxy,0.35,0.35))
cx2.plot(x1,y1,c='black')
cx2.plot(x2,y2,c='black')
cx2.plot(x3,y3,c='black')
cx2.add_patch(mpl.patches.Circle((-0.10,-0.10), radius=0.07, color=cc[0],ec='black'))
cx2.add_patch(mpl.patches.Circle((1.10,-0.10), radius=0.07, color=cc[1],ec='black'))
cx2.add_patch(mpl.patches.Circle((0.5,1.00), radius=0.07, color=cc[2],ec='black'))

#low region1
cx2.fill_between(xl1,yl1,0,color='lightgray')

#low region2
cx2.fill_between(xl2,yl2,0,color='lightgray')

#low region3
cx2.fill_between(xl3,yl3,0,color='lightgray')

#high region1
cx2.fill_between(xh1,yh1u,yh1d,color='lightgray')

#high region1
cx2.fill_between(xh2,yh2u,yh2d,color='lightgray')#,label='Artificial\n-selection\n-dominant')


cx2.scatter(i1x,i1y,color='black',label='Initial',fc='white',zorder=20)
cx2.scatter(i2x,i2y,color='black',fc='white',zorder=20)
cx2.scatter(i3x,i3y,color='black',fc='white',zorder=20)
fx=0.5
fy=np.sqrt(3)*0.2
cx2.scatter(fx,fy,color='black',label='Target')
cx2.annotate('',xy=(cv2xr,cv2yr),xytext=(i1x,i1y),arrowprops=dict(arrowstyle='->',color='red'))
cx2.annotate('',xy=(cv2xr,cv2yr),xytext=(i2x,i2y),arrowprops=dict(arrowstyle='->',color='red'))
cx2.annotate('',xy=(cv2xr,cv2yr),xytext=(i3x,i3y),arrowprops=dict(arrowstyle='->',color='red'))
cx2.axis('scaled')
cx2.axis('off')




cx3=plt.axes((cxx+0.6,cxy,0.35,0.35))
cx3.add_patch(mpl.patches.Circle((-0.10,-0.10), radius=0.07, color=cc[0],ec='black'))
cx3.add_patch(mpl.patches.Circle((1.10,-0.10), radius=0.07, color=cc[1],ec='black'))
cx3.add_patch(mpl.patches.Circle((0.5,1.00), radius=0.07, color=cc[2],ec='black'))

cx3.plot(x1,y1,c='black')
cx3.plot(x2,y2,c='black')
cx3.plot(x3,y3,c='black')


cx3.fill_between(xl1,yl1,0,color='lightgray')
cx3.fill_between(xl2,yl2,0,color='lightgray')
cx3.fill_between(xl3,yl3,0,color='lightgray')
cx3.fill_between(xh1,yh1u,yh1d,color='lightgray')
cx3.fill_between(xh2,yh2u,yh2d,color='lightgray',label='AS-dominant')



cx3.scatter(i1x,i1y,color='black',label='Initial 1',fc='white',zorder=20)
cx3.scatter(i2x,i2y,color='black',fc='white',zorder=20)
cx3.scatter(i3x,i3y,color='black',fc='white',zorder=20)
fx=0.85
fy=0.15
cx3.scatter(fx,fy,color='black',label='Target')

cx3.annotate('',xy=(cv2xr,cv2yr),xytext=(i1x,i1y),arrowprops=dict(arrowstyle='->',color='red'))
cx3.annotate('',xy=(cv1xr,cv1yr),xytext=(i2x,i2y),arrowprops=dict(arrowstyle='->',color='red'))
cx3.annotate('',xy=(fx,fy),xytext=(i3x,i3y),arrowprops=dict(arrowstyle='->',color='green'))
cx3.annotate('Success',xy=(0.4,-0.1),xytext=(0.6,-0.1),xycoords='axes fraction',arrowprops=dict(arrowstyle='<-',color='green'),va='center')
cx3.annotate('Fail',xy=(0.4,-0.22),xytext=(0.6,-0.22),xycoords='axes fraction',arrowprops=dict(arrowstyle='<-',color='red'),va='center')

cx3.axis('scaled')
cx3.axis('off')

formatter='svg' #or png
plt.savefig('figures/Fig6.'+formatter,dpi=300,bbox_inches='tight',format=formatter)

plt.show()
