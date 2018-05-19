import sys
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

Nx       = int(sys.argv[1]);
dt       = int(sys.argv[2]);
saveIter = int(sys.argv[3]);
usemixer = int(sys.argv[4]);

Ny = int(1.5*Nx);
h = 1.0/Ny;
x = np.linspace(-h/2,(1/1.5)+h/2,num=Nx+2)
y = np.linspace(1+h/2,-h/2,num=Ny+2)
XT, YT = np.meshgrid(x,y)
x = np.linspace(0,(1/1.5),num=Nx+1)
y = np.linspace(1,0,num=Ny+1)
Xv, Yv = np.meshgrid(x,y)
Tcolors = np.linspace(-0.005,0.005,num=101)
vcolors = np.linspace(0,0.05,num=101)
wcolors = np.linspace(-2,2,num=101)

cmap = plt.cm.jet
N = cmap.N
cmap.set_under(cmap(1))
cmap.set_over(cmap(N-1))

diag = open('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx,dt,usemixer),'rb');
b = os.path.getsize('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx,dt,usemixer));
diag = np.fromfile(diag,dtype=np.float64,count= int(b/8));
diag = diag.reshape((int(b/48),6));
tUH = np.linspace(1/dt,b/48/dt,num=int(b/48))
f, axarr = plt.subplots(2,3,figsize=(10,8))
plt.subplots_adjust(wspace=0.3,hspace=0.3)
axarr[0,0].plot(tUH,diag[:,0])
axarr[0,0].grid(True)
axarr[0,0].set_title('Average temperature')
axarr[0,0].set_xlabel(r'$tU/H$')
axarr[0,0].set_ylabel(r'$\frac{\langle T\rangle(t)-T_0}{\Delta T}$')
axarr[0,1].plot(tUH,diag[:,1])
axarr[0,1].grid(True)
axarr[0,1].set_title('Mixer temperature')
axarr[0,1].set_xlabel(r'$tU/H$')
axarr[0,1].set_ylabel(r'$\frac{\langle T\rangle|_{cyl}(t)-T_0}{\Delta T}$')
axarr[0,2].plot(tUH,diag[:,4])
axarr[0,2].grid(True)
axarr[0,2].set_title(r'Mesh Reynolds number $Re_h$')
axarr[0,2].set_xlabel(r'$tU/H$')
axarr[0,2].set_ylabel(r'$Re_h = \frac{(|u|+|v|)h}{\nu}$')
axarr[1,0].plot(tUH,diag[:,2])
axarr[1,0].grid(True)
axarr[1,0].set_title('RMS temperature')
axarr[1,0].set_xlabel(r'$tU/H$')
axarr[1,0].set_ylabel(r'$\frac{T_{rms}(t)-T_0}{\Delta T}$')
axarr[1,1].plot(tUH,diag[:,3])
axarr[1,1].grid(True)
axarr[1,1].set_title('Average heat flux')
axarr[1,1].set_xlabel(r'$tU/H$')
axarr[1,1].set_ylabel(r'$\frac{q_e(t)}{q_w}$')
axarr[1,2].plot(tUH,diag[:,5])
axarr[1,2].grid(True)
axarr[1,2].set_title(r'Mesh Reynolds number $Re_{h\omega}$')
axarr[1,2].set_xlabel(r'$tU/H$')
axarr[1,2].set_ylabel(r'$Re_{h\omega} = \frac{|\omega|h^2}{\nu}$')
plt.subplots_adjust(top=0.8)
plt.savefig('results/diagnostics_Nx%d_dt%d_mixing%d.eps' % (Nx,dt,usemixer))
plt.show()

def framesT():
    return len(glob.glob('data/Nx%d_dt%d_mixing%d/T_iter*.bin' % (Nx,dt,usemixer)))

def framesv():
    return len(glob.glob('data/Nx%d_dt%d_mixing%d/v_iter*.bin' % (Nx,dt,usemixer)))

def framesw():
    return len(glob.glob('data/Nx%d_dt%d_mixing%d/w_iter*.bin' % (Nx,dt,usemixer)))

def initT():
	fd = open('data/Nx%d_dt%d_mixing%d/T_iter%d.bin' % (Nx,dt,usemixer,0),'rb');
	T = np.fromfile(fd,dtype=np.float64, count = (Nx+2)*(Ny+2));
	T = T.reshape((Ny+2,Nx+2));
	CS = plt.contourf(XT,YT,T,Tcolors,cmap=cmap,extend="both")
	
	cbar = plt.colorbar(CS)
	plotmixer(usemixer,0)
	return CS

def initv():
    fd = open('data/Nx%d_dt%d_mixing%d/v_iter%d.bin' % (Nx,dt,usemixer,0),'rb');
    v = np.fromfile(fd,dtype=np.float64, count = (Nx+1)*(Ny+1));
    v = v.reshape((Ny+1,Nx+1));

    CS = plt.contourf(Xv,Yv,v,vcolors,cmap=cmap,extend="both")
    cbar = plt.colorbar(CS)
    plotmixer(usemixer,0)
    return CS

def initw():
    fd = open('data/Nx%d_dt%d_mixing%d/w_iter%d.bin' % (Nx,dt,usemixer,0),'rb');
    w = np.fromfile(fd,dtype=np.float64, count = (Nx+1)*(Ny+1));
    w = w.reshape((Ny+1,Nx+1));

    CS = plt.contourf(Xv,Yv,w,wcolors,cmap=cmap,extend="both")
    cbar = plt.colorbar(CS)
    plotmixer(usemixer,0)

    return CS

def plotmixer(um,i):
    if (um == 1):
        theta = np.linspace(0,2*np.pi,1000)
        omega = 0.1
        r = 0.2*np.cos(3*(theta - i*saveIter*omega/dt))
        x = r*np.cos(theta) + 1/3
        y = r*np.sin(theta) + 1/3
        ax.fill_between(x,0,y,facecolor='black')
        x = 0.04*np.cos(theta) + 1/3
        y = 0.04*np.sin(theta) + 1/3
        ax.fill_between(x,0,y,facecolor='black')


def animateT(i,CS,mf):
    print('T: %d/%d (%.1f%%)\r' % (i*saveIter, mf*saveIter, i/mf*100),end="",flush=True)
    ax.clear()
    fd = open('data/Nx%d_dt%d_mixing%d/T_iter%d.bin' % (Nx,dt,usemixer,i*saveIter),'rb');
    T = np.fromfile(fd,dtype=np.float64, count = (Nx+2)*(Ny+2));
    try:
        T = T.reshape((Ny+2,Nx+2));
    except:
        return CS
    CS = plt.contourf(XT,YT,T,Tcolors,cmap=cmap,extend="both")
    plotmixer(usemixer,i)
    plt.title(r'$\frac{T-T_0}{\Delta T}$ at $\frac{tU}{H} = %.2f$' % (i*saveIter/dt))
    return CS

def animatev(i,CS,mf):
    print('v: %d/%d (%.1f%%)\r' % (i*saveIter, mf*saveIter, i/mf*100),end="",flush=True)
    ax.clear()
    fd = open('data/Nx%d_dt%d_mixing%d/v_iter%d.bin' % (Nx,dt,usemixer,i*saveIter),'rb');
    v = np.fromfile(fd,dtype=np.float64, count = (Nx+1)*(Ny+1));
    try:
        v = v.reshape((Ny+1,Nx+1));
    except:
        return CS
    CS = plt.contourf(Xv,Yv,v,vcolors,cmap=cmap,extend="both")
    plotmixer(usemixer,i)
    plt.title(r'$\frac{|v|}{U}$ at $\frac{tU}{H} = %.2f$' % (i*saveIter/dt))
    return CS

def animatew(i,CS,mf):
    print('w: %d/%d (%.1f%%)\r' % (i*saveIter, mf*saveIter, i/mf*100),end="",flush=True)
    ax.clear()
    fd = open('data/Nx%d_dt%d_mixing%d/w_iter%d.bin' % (Nx,dt,usemixer,i*saveIter),'rb');
    w = np.fromfile(fd,dtype=np.float64, count = (Nx+1)*(Ny+1));
    try:
        w = w.reshape((Ny+1,Nx+1));
    except:
        return CS
    CS = plt.contourf(Xv,Yv,w,wcolors,cmap=cmap,extend="both")
    plotmixer(usemixer,i)
    plt.title(r'$\frac{\omega H}{U}$ at $\frac{tU}{H} = %.2f$' % (i*saveIter/dt))
    return CS


fig = plt.figure()
ax = plt.axes(xlim=(0,1/1.5),ylim=(0,1))
CS = initT()
maxframe = framesT()
anim = animation.FuncAnimation(fig,animateT,interval=1000/25,fargs=(CS,maxframe,),frames=maxframe)
anim.save('results/T_Nx%d_dt%d_mixing%d.mp4' % (Nx,dt,usemixer))
plt.close(fig)
print('T: done                    ')

fig = plt.figure()
ax = plt.axes(xlim=(0,1/1.5),ylim=(0,1))
CS = initv()
maxframe = framesv()
anim = animation.FuncAnimation(fig,animatev,interval=1000/25,fargs=(CS,maxframe,),frames=maxframe)
anim.save('results/v_Nx%d_dt%d_mixing%d.mp4' % (Nx,dt,usemixer))
plt.close(fig)
print('v: done                      ')

fig = plt.figure()
ax = plt.axes(xlim=(0,1/1.5),ylim=(0,1))
CS = initw()
maxframe = framesw()
anim = animation.FuncAnimation(fig,animatew,interval=1000/25,fargs=(CS,maxframe,),frames=maxframe)
anim.save('results/w_Nx%d_dt%d_mixing%d.mp4' % (Nx,dt,usemixer))
plt.close(fig)
print('w: done                        ')

