import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpltern


###############################3
# Draw three strain system
###############################


plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'DejaVu Sans:italic:bold'
plt.rcParams['mathtext.bf'] = 'DejaVu Sans:italic:bold'
plt.rcParams['font.family'] = 'Arial'
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

#parameter prepare
mu=1e-4
r=0.5
s=3e-2
N0=1000
mbar=100
ncomm=10
g=ncomm
rhat=0
nensemble=30
ncycle=30
tcycle=np.log(ncomm+1)/r

#######################
# Define accessible region
#######################

data=np.loadtxt("N0%d_r%s_s%s_mu%s_ext_means_fine.dat"%(N0,r,s,mu))
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


color='#d8b365ff'
oldgrey=mpl.cm.get_cmap('Greys')
newgrey=mpl.colors.ListedColormap(['white','#d8b365ff'])#['#b4b4b4ff','#d8b365ff'])#oldgrey(np.linspace(0,0.5,3)))
#tc=ax.tripcolor(fvs,fws,fms,dscore,shading='gouraud',rasterized=True,cmap=newgrey)

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
    
def draw_background(ax):
    #draw background
    ax.tripcolor(fvs,fws,fms,dscore,shading='gouraud',rasterized=True,cmap=newgrey)
    return ax

def Draw_triangle_from_data(ax,mvs,mvhats,tgtcol='black',**kwargs):
    ###
    # Draw composition trajectory in triangle plot with accesible region
    # mvs =[ (m0,v0,trjcol)...]
    ###

    for (mv,mvhat) in itt.product(mvs,mvhats):
    
        #Initia/Target point
        m0,v0,trjcol=mv
        mhat,vhat=mvhat
        fm=m0/N0
        fv=v0/N0
        fw=1-fm-fv
        #ax.scatter(fv,fw,fm,marker='o',s=20,ec='black',fc='whitle',zorder=10,label='Initial')
        ax.scatter(fv,fw,fm,marker='o',s=20,ec='black',fc='white',zorder=10,label='Initial')
        ax.scatter(vhat,1-vhat-mhat,mhat,marker='^',s=30,c=tgtcol,zorder=10,label='Target')
    
        AGS=[]
        folder="data/raw/"
        nensemble=30
        for e in range(nensemble):
            try:
                AGS.append(np.loadtxt(folder+"tri/N0%s_m0%s_v0%s_r%s_s%s_mu%s_g%s_mhat%s_vhat%s_AGS_point_%d.cycle"%(N0,m0,v0,r,s,mu,ncomm,mhat,vhat,e)))
            except:
                AGS.append(np.loadtxt(folder+"AGS_PD_sto_N0%s_mbar%s_vbar%s_r%s_s%s_mu%s_ncomm%s_mhat%s_vhat%s_ncycle%d/%d.cycle"%(N0,m0,v0,r,s,mu,ncomm,mhat,vhat,ncycle,e)))
        AGS=np.array(AGS) # [ensemble, timestamp, (T,selected index,w----,m----,v----)]
        
        T=AGS[0,:,0]
        sel_inds=AGS[:,:,1].astype(int) #selected index - [ensemble, timestamp]
        c1_sel=np.zeros((nensemble,len(T),ncomm)) #c avg datas, [ensemble, timestamp]
        c2_sel=np.zeros((nensemble,len(T),ncomm)) #c avg datas, [ensemble, timestamp]
        w_cho=np.zeros((nensemble,len(T))) #w datas, [ensemble, timestamp]
        m_cho=np.zeros((nensemble,len(T))) #m datas, [ensemble, timestamp]
        v_cho=np.zeros((nensemble,len(T))) #m datas, [ensemble, timestamp]
        #ensemble for initial
        for i in range(nensemble): #timestamp
            for j in range(len(T)):
                #c_sel[i,j,:]=cf(AGS[i,j,ncomm+2:ncomm*2+2],AGS[i,j,ncomm*3+2:])
                #w_cho[i,j]=AGS[i,j,sel_inds[i,j]+2+ncomm]
                #m_cho[i,j]=AGS[i,j,sel_inds[i,j]+2+ncomm*3]
                #v_cho[i,j]=AGS[i,j,sel_inds[i,j]+2+ncomm*5]
                w_cho[i,j]=AGS[i,j,sel_inds[i,j]+2]
                m_cho[i,j]=AGS[i,j,sel_inds[i,j]+2+ncomm*2]
                v_cho[i,j]=AGS[i,j,sel_inds[i,j]+2+ncomm*4]
        
        #Artificial group selection with the selected
        c0_cho=cf0(w_cho,m_cho,v_cho)
        c1_cho=cf1(w_cho,m_cho,v_cho)
        c2_cho=cf2(w_cho,m_cho,v_cho)
        cho0avg=np.mean(c0_cho,axis=0)    
        cho1avg=np.mean(c1_cho,axis=0)    
        cho2avg=np.mean(c2_cho,axis=0)
        #lines=ax.plot(cho2avg,cho0avg,cho1avg)    
        #print(lines)
        alphas=np.linspace(0.1,1,len(cho2avg[::1]))
        #for i in range(len(alphas)):
        #    print(lines[i])
        #    lines[i].set(alpha=alphas[i])
        
        ax.scatter(cho2avg[::1],cho0avg[::1],cho1avg[::1],marker='o',alpha=alphas,color=trjcol,s=10)    
        #[ax.plot(cho2avg[i:i+2],cho0avg[i:i+2],cho1avg[i:i+2],c=trjcol,alpha=alphas[i]) for i in range(len(cho2avg)-1)]

    return ax

def Draw_triangle_from_data_single(ax,mvs,mvhats,lidx,tgtcol='black'):
    ###
    # Draw composition trajectory in triangle plot with accesible region
    # mvs =[ (m0,v0,trjcol)...]
    # lidx=[index of lines]
    ###

    for (mv,mvhat) in itt.product(mvs,mvhats):
    
        #Initia/Target point
        m0,v0,trjcol=mv
        mhat,vhat=mvhat
        fm=m0/N0
        fv=v0/N0
        fw=1-fm-fv
        #ax.scatter(fv,fw,fm,marker='o',s=20,ec='black',fc='whitle',zorder=10,label='Initial')
        ax.scatter(fv,fw,fm,marker='o',s=20,c='black',zorder=10,label='Initial')
        ax.scatter(vhat,1-vhat-mhat,mhat,marker='^',s=30,c=tgtcol,zorder=10,label='Target')
    
        AGS=[]
        folder="data/raw/"
        nensemble=30
        for e in range(nensemble):
            #AGS.append(np.loadtxt(folder+"N0%s_m0%s_v0%s_r%s_s%s_mu%s_g%s_mhat%s_vhat%s_AGS_point_%d.cycle"%(N0,m0,v0,r,s,mu,ncomm,mhat,vhat,e)))
            AGS.append(np.loadtxt(folder+"AGS_PD_sto_N0%s_mbar%s_vbar%s_r%s_s%s_mu%s_ncomm%s_mhat%s_vhat%s_ncycle%d%d.cycle"%(N0,m0,v0,r,s,mu,ncomm,mhat,vhat,ncycle,e)))
        AGS=np.array(AGS) # [ensemble, timestamp, (T,selected index,w----,m----,v----)]
        
        T=AGS[0,:,0]
        sel_inds=AGS[:,:,1].astype(int) #selected index - [ensemble, timestamp]
        c1_sel=np.zeros((nensemble,len(T),ncomm)) #c avg datas, [ensemble, timestamp]
        c2_sel=np.zeros((nensemble,len(T),ncomm)) #c avg datas, [ensemble, timestamp]
        w_cho=np.zeros((nensemble,len(T))) #w datas, [ensemble, timestamp]
        m_cho=np.zeros((nensemble,len(T))) #m datas, [ensemble, timestamp]
        v_cho=np.zeros((nensemble,len(T))) #m datas, [ensemble, timestamp]
        #ensemble for initial
        for i in range(nensemble): #timestamp
            for j in range(len(T)):
                #c_sel[i,j,:]=cf(AGS[i,j,ncomm+2:ncomm*2+2],AGS[i,j,ncomm*3+2:])
                #w_cho[i,j]=AGS[i,j,sel_inds[i,j]+2+ncomm]
                #m_cho[i,j]=AGS[i,j,sel_inds[i,j]+2+ncomm*3]
                #v_cho[i,j]=AGS[i,j,sel_inds[i,j]+2+ncomm*5]
                w_cho[i,j]=AGS[i,j,sel_inds[i,j]+2]
                m_cho[i,j]=AGS[i,j,sel_inds[i,j]+2+ncomm*2]
                v_cho[i,j]=AGS[i,j,sel_inds[i,j]+2+ncomm*4]
    
        #Artificial group selection with the selected
        c0_cho=cf0(w_cho,m_cho,v_cho)
        c1_cho=cf1(w_cho,m_cho,v_cho)
        c2_cho=cf2(w_cho,m_cho,v_cho)
        alphas=np.linspace(0.1,1,len(c2_cho[0,::1]))
        for i in lidx:    
            ax.scatter(c2_cho[i,::1],c0_cho[i,::1],c1_cho[i,::1],marker='o',color='green',s=10,alpha=alphas)    

    return ax
succ='black'#'green'
fail='black'#'red'

cxx=mxx+0.35#0.55
cxy=mxy+0.1#0.37
cx=plt.axes((cxx,cxy,0.30,0.30),projection='ternary')
cx.annotate('b',xy=(-0.0,1.1),weight='bold',xycoords='axes fraction')
cx.annotate(' Simulation of composition trajectories according to initial/target points',xy=(0.05,1.1),xycoords='axes fraction')
cx=draw_background(cx)
cx=Draw_triangle_from_data(cx,[(50,50,succ),(150,5,succ)],[(0.02,0.9)],tgtcol='#d62728')
cx.taxis.set_ticks([])
cx.raxis.set_ticks([])
cx.laxis.set_ticks([])

cx2=plt.axes((cxx+0.35,cxy,0.30,0.30),projection='ternary')
cx2=draw_background(cx2)
cx2=Draw_triangle_from_data(cx2,[(50,50,fail),(150,5,fail),(75,900,fail)],[(0.33,0.33)],tgtcol='#2ca02c')
cx2.taxis.set_ticks([])
cx2.raxis.set_ticks([])
cx2.laxis.set_ticks([])

cx3=plt.axes((cxx+0.70,cxy,0.30,0.30),projection='ternary')
cx3=draw_background(cx3)
cx3=Draw_triangle_from_data(cx3,[(50,50,fail),(150,5,succ),(75,900,fail)],[(0.75,0.05)],tgtcol='#ff00ff')
#cx3=Draw_triangle_from_data_single(cx3,[(150,5,succ)],[(0.5,0.1)],range(nensemble))
cx3.taxis.set_ticks([])
cx3.raxis.set_ticks([])
cx3.laxis.set_ticks([])

#Draw effective rest-of-system line
d=0.01
ffs=np.arange(0,1+0.5*d,d)

for a in np.arange(0.2,1,0.2):
    ss=(1-a)-(1-a)*ffs
    fs=1-ffs-ss

    cx2.plot(ffs,ss,fs,alpha=0.4,color='black',ls='--')




formatter='svg'
#plt.tight_layout()
#plt.savefig('figures/Fig5_v3_fine.'+formatter,bbox_inches='tight',dpi=300,format=formatter)
plt.show()

