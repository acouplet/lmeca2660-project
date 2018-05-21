import sys
import os
import glob
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


saveIter = 20


Nx = [150,200,200,200,300,400]
dt = [100,100,80 ,60 ,100,100]
style = ['r','b','g','r--','b--','g--']

fig = plt.figure(figsize=(8,6))
ax = plt.axes()
axins = zoomed_inset_axes(ax, 2, loc=4) # zoom-factor: 2.5, location: upper-left

for i in range(0,6):
	usemixer = 0;
	diag = open('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx[i],dt[i],usemixer),'rb');
	b = os.path.getsize('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx[i],dt[i],usemixer));
	b = int(np.floor(b/48)*48)
	
	diag = np.fromfile(diag,dtype=np.float64,count= int(b/8));
	diag = diag.reshape((int(b/48),6));
	tUH = np.linspace(1/dt[i],b/48/dt[i],num=int(b/48))
	
	ax.plot(tUH,diag[:,0],style[i],label=(r'$N_x = %d$ and $\Delta t = %.4lf$' % (Nx[i],1/dt[i])))
	axins.plot(tUH,diag[:,0],style[i],alpha=0.5)
	
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::1],labels[::1])
ax.grid(True)
ax.set_title('Evolution of the average temperature without mixer')
ax.set_xlabel(r'$tU/H$')
ax.set_ylabel(r'$\frac{\langle T\rangle(t)-T_0}{\Delta T}$')
axins.set_xlim(100,600)
axins.set_ylim(0.0003,0.001)	
axins.grid(True)
plt.yticks(visible=True)
axins.xaxis.tick_top()
plt.xticks(visible=True)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
plt.savefig('results/all_Tavg_mixing%d.pdf' % (usemixer))


fig = plt.figure(figsize=(8,6))
ax = plt.axes()
axins = zoomed_inset_axes(ax, 2, loc=4) # zoom-factor: 2.5, location: upper-left

for i in range(0,6):
	usemixer = 1;
	diag = open('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx[i],dt[i],usemixer),'rb');
	b = os.path.getsize('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx[i],dt[i],usemixer));
	b = int(np.floor(b/48)*48)
	
	diag = np.fromfile(diag,dtype=np.float64,count= int(b/8));
	diag = diag.reshape((int(b/48),6));
	tUH = np.linspace(1/dt[i],b/48/dt[i],num=int(b/48))
	
	ax.plot(tUH,diag[:,0],style[i],label=(r'$N_x = %d$ and $\Delta t = %.4lf$' % (Nx[i],1/dt[i])))
	axins.plot(tUH,diag[:,0],style[i],alpha=0.5)
	
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::1],labels[::1])
ax.grid(True)
ax.set_title('Evolution of the average temperature with mixer')
ax.set_xlabel(r'$tU/H$')
ax.set_ylabel(r'$\frac{\langle T\rangle(t)-T_0}{\Delta T}$')
axins.set_xlim(100,220)
axins.set_ylim(0.0002,0.0005)	
axins.grid(True)
plt.yticks(visible=True)
axins.xaxis.tick_top()
plt.xticks(visible=True)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
plt.savefig('results/all_Tavg_mixing%d.pdf' % (usemixer))

fig = plt.figure(figsize=(8,6))
ax = plt.axes()

for i in range(0,6):
	usemixer = 1;
	diag = open('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx[i],dt[i],usemixer),'rb');
	b = os.path.getsize('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx[i],dt[i],usemixer));
	b = int(np.floor(b/48)*48)
	
	diag = np.fromfile(diag,dtype=np.float64,count= int(b/8));
	diag = diag.reshape((int(b/48),6));
	tUH = np.linspace(1/dt[i],b/48/dt[i],num=int(b/48))
	
	ax.plot(tUH,diag[:,1],style[i],label=(r'$N_x = %d$ and $\Delta t = %.4lf$' % (Nx[i],1/dt[i])))
	
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::1],labels[::1])
ax.grid(True)
ax.set_title("Evolution of the mixer's average temperature")
ax.set_xlabel(r'$tU/H$')
ax.set_ylabel(r'$\frac{\langle T\rangle(t)-T_0}{\Delta T}$')
plt.savefig('results/all_TMavg_mixing%d.pdf' % (usemixer))


##############
## REYNOLDS ##
##############

Nx = [200,200,200,300,400]
dt = [100,80 ,60 ,100,100]
style = ['b','g','r--','b--','g--']

fig = plt.figure(figsize=(8,6))
ax = plt.axes()

for i in range(0,len(Nx)):
	usemixer = 0;
	diag = open('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx[i],dt[i],usemixer),'rb');
	b = os.path.getsize('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx[i],dt[i],usemixer));
	b = int(np.floor(b/48)*48)
	
	diag = np.fromfile(diag,dtype=np.float64,count= int(b/8));
	diag = diag.reshape((int(b/48),6));
	tUH = np.linspace(1/dt[i],b/48/dt[i],num=int(b/48))
	
	ax.plot(tUH,diag[:,4],style[i],label=(r'$N_x = %d$ and $\Delta t = %.4lf$' % (Nx[i],1/dt[i])))
	
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::1],labels[::1])
ax.grid(True)
ax.set_title(r'Mesh Reynolds number $Re_h$ without mixer')
ax.set_xlabel(r'$tU/H$')
ax.set_ylabel(r'$Re_h = \frac{(|u|+|v|)h}{\nu}$')
plt.savefig('results/all_Reh_mixing%d.pdf' % (usemixer))

fig = plt.figure(figsize=(8,6))
ax = plt.axes()

for i in range(0,len(Nx)):
	usemixer = 1;
	diag = open('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx[i],dt[i],usemixer),'rb');
	b = os.path.getsize('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx[i],dt[i],usemixer));
	b = int(np.floor(b/48)*48)
	
	diag = np.fromfile(diag,dtype=np.float64,count= int(b/8));
	diag = diag.reshape((int(b/48),6));
	tUH = np.linspace(1/dt[i],b/48/dt[i],num=int(b/48))
	
	ax.plot(tUH,diag[:,4],style[i],label=(r'$N_x = %d$ and $\Delta t = %.4lf$' % (Nx[i],1/dt[i])))
	
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::1],labels[::1])
ax.grid(True)
ax.set_title(r'Mesh Reynolds number $Re_h$ with mixer')
ax.set_xlabel(r'$tU/H$')
ax.set_ylabel(r'$Re_h = \frac{(|u|+|v|)h}{\nu}$')
plt.savefig('results/all_Reh_mixing%d.pdf' % (usemixer))

fig = plt.figure(figsize=(8,6))
ax = plt.axes()

for i in range(0,len(Nx)):
	usemixer = 0;
	diag = open('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx[i],dt[i],usemixer),'rb');
	b = os.path.getsize('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx[i],dt[i],usemixer));
	b = int(np.floor(b/48)*48)
	
	diag = np.fromfile(diag,dtype=np.float64,count= int(b/8));
	diag = diag.reshape((int(b/48),6));
	tUH = np.linspace(1/dt[i],b/48/dt[i],num=int(b/48))
	
	ax.plot(tUH,diag[:,5],style[i],label=(r'$N_x = %d$ and $\Delta t = %.4lf$' % (Nx[i],1/dt[i])))
	
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::1],labels[::1])
ax.grid(True)
ax.set_title(r'Mesh Reynolds number $Re_{h,\omega}$ without mixer')
ax.set_xlabel(r'$tU/H$')
ax.set_ylabel(r'$Re_{h\omega} = \frac{|\omega|h^2}{\nu}$')
plt.savefig('results/all_Rehw_mixing%d.pdf' % (usemixer))

fig = plt.figure(figsize=(8,6))
ax = plt.axes()

for i in range(0,len(Nx)):
	usemixer = 1;
	diag = open('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx[i],dt[i],usemixer),'rb');
	b = os.path.getsize('data/Nx%d_dt%d_mixing%d/diagnostics.bin' % (Nx[i],dt[i],usemixer));
	b = int(np.floor(b/48)*48)
	
	diag = np.fromfile(diag,dtype=np.float64,count= int(b/8));
	diag = diag.reshape((int(b/48),6));
	tUH = np.linspace(1/dt[i],b/48/dt[i],num=int(b/48))
	
	ax.plot(tUH,diag[:,5],style[i],label=(r'$N_x = %d$ and $\Delta t = %.4lf$' % (Nx[i],1/dt[i])))
	
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::1],labels[::1])
ax.grid(True)
ax.set_title(r'Mesh Reynolds number $Re_{h,\omega}$ with mixer')
ax.set_xlabel(r'$tU/H$')
ax.set_ylabel(r'$Re_{h\omega} = \frac{|\omega|h^2}{\nu}$')
plt.savefig('results/all_Rehw_mixing%d.pdf' % (usemixer))

plt.show()
	
