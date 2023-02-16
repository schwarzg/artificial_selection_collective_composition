import matplotlib.pyplot as plt
import numpy as np
#Two strain case
ax=plt.axes((0.1,0.7,0.4,0.3))


#ax.hlines(0,0,1,colors='black',lw=3)
ax.annotate('',xy=(1,0),xytext=(0,0),arrowprops=dict(arrowstyle='->'))
#ax.fill_between([0.3,0.7],0.02,0,color='lightgray')
ax.fill_between([0,0.3],0.02,0,color='green',alpha=0.5)
ax.fill_between([0.7,1],0.02,0,color='red',alpha=0.5)
#ax.annotate('Fail',xy=(0.45,0.005))
ax.annotate(r'$0$',xy=(0,-0.02))
ax.annotate(r'$1$',xy=(0.95,-0.02))
ax.annotate('Mutant\nfrequency',xy=(1.05,-0.005))
ax.annotate('(a) Two-strain',xy=(-0.2,0.95),xycoords='axes fraction')

ax.set_ylim(ymin=-0.05,ymax=0.05)

ax.scatter(0.25,0.01,c='black',marker='s')
ax.annotate('Target',xy=(0.25,0.015),xytext=(0.2,0.03),arrowprops=dict(arrowstyle='->'))
ax.scatter(0.9,0.01,c='black')
ax.annotate('Initial',xy=(0.9,0.015),xytext=(0.8,0.03),arrowprops=dict(arrowstyle='->'))
ax.scatter(0.05,0.01,c='black',marker='v')
ax.annotate('Initial',xy=(0.05,0.015),xytext=(-0.05,0.03),arrowprops=dict(arrowstyle='->'))
ax.annotate('',xy=(0.25,0.01),xytext=(0.05,0.01),arrowprops=dict(arrowstyle='->'))
ax.annotate('',xy=(0.25,0.01),xytext=(0.9,0.01),arrowprops=dict(arrowstyle='->'))
#ax.plot([0.25,0.9],[0.01,0.01],ls=':',color='black')
ax.scatter([0.5],[0.01],marker='x',color='red')
ax.annotate('Fail',xy=(0.5,0.01),xytext=(0.45,-0.03),arrowprops=dict(arrowstyle='->'))
ax.scatter([0.15],[0.01],facecolor='none',edgecolor='blue')
ax.annotate('Success',xy=(0.15,0.01),xytext=(0.05,-0.03),arrowprops=dict(arrowstyle='->'))

ax.vlines([0.3,0.7],0.005,-0.005,colors='black')
ax.annotate(r'$f^L$',xy=(0.28,-0.02))
ax.annotate(r'$f^U$',xy=(0.68,-0.02))

ax.legend(frameon=False,loc=(0,-0.1))
ax.axis('off')


#Three strain case
bx=plt.axes((0.1,0.2,0.4,0.4))
bx.annotate('(b) Three-strain',xy=(-0.2,1.1),xycoords='axes fraction')

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
cv1y=curve1(cv1x,a1,b1)

cv2yl=0.71
cv2xl=cv2yl/np.sqrt(3)
cv2yr=0.7
cv2xr=1-cv2yr/np.sqrt(3)

def curve2(x,a,b,c):
	return a*(x-c)**2+b

c=0.51
a2=(cv2yl-cv2yr)/((cv2xl-c)**2-(cv2xr-c)**2)
b2=(-(cv2xr-c)**2*cv2yl+(cv2xl-c)**2*cv2yr)/((cv2xl-c)**2-(cv2xr-c)**2)

cv2x=np.arange(cv2xl,cv2xr+0.005,0.01)
cv2y=curve2(cv2x,a2,b2,c)

bx.plot(x1,y1,c='black')
bx.plot(x2,y2,c='black')
bx.plot(x3,y3,c='black')
bx.annotate('Wild type\nfrequency',xy=(-0.1,-0.15),xycoords='axes fraction')
bx.annotate('Single mutant\nfrequency',xy=(0.83,-0.15),xycoords='axes fraction')
bx.annotate('Double mutant\nfrequency',xy=(0.3,1.02),xycoords='axes fraction')


#divide 3regions

#region 1
#x1=np.linspace(cv1xl,cv2xl,30)
#y1u=x1*np.sqrt(3)
#y1d=curve1(x1,a1,b1)
#bx.fill_between(x1,y1u,y1d,color='lightgray')

#region2
#x2=np.linspace(cv2xl,cv2xr,30)
#y2u=curve2(x2,a2,b2,c)
#y2d=curve1(x2,a1,b1)
#bx.fill_between(x2,y2u,y2d,color='lightgray')
#bx.annotate('Fail',xy=(0.45,0.4))

#region2
#x3=np.linspace(cv2xr,cv1xr,30)
#y3u=-np.sqrt(3)*(x3-1)
#y3d=curve1(x3,a1,b1)
#bx.fill_between(x3,y3u,y3d,color='lightgray')

#low region1
xl1=np.linspace(0,cv1xl)
yl1=np.sqrt(3)*xl1
bx.fill_between(xl1,yl1,0,color='lightgray')

#low region2
xl2=np.linspace(cv1xl,cv1xr)
#yl2=curve1(xl2,a1,b1)
yl2=curve1p(xl2,cv1xl,cv1yl,cv1xr,cv1yr)
bx.fill_between(xl2,yl2,0,color='lightgray')

#low region3
xl3=np.linspace(cv1xr,1)
yl3=-np.sqrt(3)*(xl3-1)
bx.fill_between(xl3,yl3,0,color='lightgray')

#high region1
xh1=np.linspace(cv2xl,0.5)
yh1u=np.sqrt(3)*xh1
yh1d=curve2(xh1,a2,b2,c)
bx.fill_between(xh1,yh1u,yh1d,color='lightgray')

#high region1
xh2=np.linspace(0.5,cv2xr)
yh2u=-np.sqrt(3)*(xh2-1)
yh2d=curve2(xh2,a2,b2,c)
bx.fill_between(xh2,yh2u,yh2d,color='lightgray',label='Artificial\n-selection\n-dominant')


#point initial and target
i1x=0.5
i1y=0.75
bx.scatter(i1x,i1y,color='black',label='Initial 1')
bx.annotate('Initial',xy=(i1x,i1y),xytext=(i1x+0.1,i1y),arrowprops=dict(arrowstyle='->'))



i2x=0.05
i2y=0.03
bx.scatter(i2x,i2y,color='black',marker='v',label='Initial 2')
bx.annotate('Initial',xy=(i2x-0.01,i2y-0.01),xytext=(i2x-0.15,i2y+0.1),arrowprops=dict(arrowstyle='->'))

i3x=0.3
i3y=0.00
bx.scatter(i3x,i3y,color='black',marker='^',label='Initial 3')
bx.annotate('Initial',xy=(i3x,i3y-0.01),xytext=(i3x,i3y-0.2),arrowprops=dict(arrowstyle='->'))

fx=0.85
fy=0.15
bx.scatter(fx,fy,marker='s',color='black',label='Target')
bx.annotate('Target',xy=(fx,fy),xytext=(fx+0.1,fy),arrowprops=dict(arrowstyle='->'))


#for x1 to f
xline=np.arange(i1x,fx,0.01)
path1=lambda x:(x-fx)*(i1y-fy)/(i1x-fx)+fy
diff=lambda x:curve2(x,a2,b2,c)-path1(x)
import scipy.optimize as opt

crossx=opt.fsolve(diff,0.5)[0]
crossy=curve2(crossx,a2,b2,c)

bx.annotate('',xy=(fx,fy),xytext=(i1x,i1y),arrowprops=dict(arrowstyle='->'))
#bx.plot([i1x,fx],[i1y,fy],c='black',ls=':')
locx1x=0.5*(fx+i1x)
locx1y=0.5*(i1y+fy)
bx.scatter([locx1x],[locx1y],marker='x',color='red')
bx.annotate('Fail',xy=(locx1x,locx1y),xytext=(locx1x+0.2,locx1y),arrowprops=dict(arrowstyle='->'))

#for x2 to f
bx.annotate('',xy=(fx,fy),xytext=(i2x,i2y),arrowprops=dict(arrowstyle='->'))
bx.scatter([0.5*(fx+i2x)],[0.5*(i2y+fy)],marker='x',color='red')
locx2x=0.5*(fx+i2x)
locx2y=0.5*(i2y+fy)
bx.annotate('Fail',xy=(locx2x,locx2y),xytext=(locx2x-0.1,locx2y+0.2),arrowprops=dict(arrowstyle='->'))

#for x3 to f
bx.annotate('',xy=(fx,fy),xytext=(i3x,i3y),arrowprops=dict(arrowstyle='->'))
bx.scatter([0.5*(fx+i3x)],[0.5*(i3y+fy)],facecolor='none',edgecolor='blue')
locx3x=0.5*(fx+i3x)
locx3y=0.5*(i3y+fy)
bx.annotate('Success',xy=(locx3x,locx3y),xytext=(locx3x,locx3y-0.25),arrowprops=dict(arrowstyle='->'))

#bx.legend(frameon=False,loc=(0.8,0.3))
bx.axis('off')

formatter='svg'
#plt.savefig('Fig6.'+formatter,dpi=300,bbox_inches='tight',format=formatter)

plt.show()
