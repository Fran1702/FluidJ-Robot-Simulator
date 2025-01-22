#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 10:37:54 2024

@author: hector.ortiz
"""
import forces_eq
import manip
import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#%%
def anim_xy(x,y):
    fig, ax = plt.subplots()
    xdata, ydata = [], []
    ln, = plt.plot([], [], 'r')
    point, = ax.plot([], [], ls="none", marker="o")
    
    def init():
        dxy = np.max((np.abs(np.max(x)-np.min(x)),np.abs(np.max(y)-np.min(y))))*1.1
        x_c = (np.max(x)+np.min(x))/2
        y_c = (np.max(y)+np.min(y))/2
        ax.set_xlim(x_c-dxy,x_c+dxy)
        ax.set_ylim(y_c-dxy,y_c+dxy)
        ln.set_data(xdata,ydata)
        point.set_data(xdata,ydata)
        return ln,point
    
    def update(k):
        i = min(k, x.size)
        #print(x[:i])
        ln.set_data(x[:i], y[:i])
        point.set_data(x[i], y[i])
        return ln,point
    
    ani = animation.FuncAnimation(fig=fig,func=update, frames=range(x.size), 
                                  init_func=init,interval=100, blit=True)
    ani.save(filename="anim.mp4", dpi =80, fps=20)
    plt.show()	

#%% ------------- Robot tested experimentally -----------------
#P_end = np.array([0,0,630])
P_end = np.array([-380,0,630]) # ANTENNA as end effector
#P_end = np.array([0,0,0])
#DROPLET_VOLUME = 0.12e-9#0.22e-9
#DROPLET_VOLUME = 0.22e-9
DROPLET_VOLUME = 0.102e-9#0.22e-9
theta_0 = 135*np.pi/180
theta_min = 60*np.pi/180
r_min = (3*DROPLET_VOLUME*np.sin(theta_0)**2/(np.pi*(2+np.cos(theta_0))*(1-np.cos(theta_0))**2))**(1/3)*1e6
r_max = (3*DROPLET_VOLUME*np.sin(theta_min)**2/(np.pi*(2+np.cos(theta_min))*(1-np.cos(theta_min))**2))**(1/3)*1e6
#r_min = (DROPLET_VOLUME*3/(2*np.pi))**(1/3)*1e6
r_min = np.round(r_min*1.0,0)
r_max = np.round(r_max*1.052,0)
print(r_min)
print(r_max)
#robot = manip.Manip(P_end=P_end, V0=DROPLET_VOLUME,R_DROPLET_MIN=480, R_DROPLET_MAX=540,
robot = manip.Manip(P_end=P_end, V0=DROPLET_VOLUME,R_DROPLET_MIN=r_min, R_DROPLET_MAX=r_max,                    
                    ZMAX=600, ZMIN=200, R_top=150)
# Load Data
robot.load_data()


 
#%% Solve points for forward kinematics
#rb_l = np.arange(robot.R_DROPLET_MIN*1e-6, robot.R_DROPLET_MAX*1e-6+1e-8, 50e-6) 
rb_l = np.linspace(robot.R_DROPLET_MIN*1e-6, robot.R_DROPLET_MAX*1e-6,8)
rb_l = np.around(rb_l,6)
# Create mesh grid
rb_arr = np.array(np.meshgrid(rb_l,rb_l,rb_l)).T.reshape(-1,3)

#rb_arr = np.array(np.meshgrid(rb_l,min(rb_l),rb_l)).T.reshape(-1,3)

print(rb_arr.shape)
# plot grid

pts = rb_arr

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(rb_arr[:,0]*1e6,rb_arr[:,1]*1e6,rb_arr[:,2]*1e6, alpha=0.8,s=2)
ax.set_xlabel('Rb1')
ax.set_ylabel('Rb2')
ax.set_zlabel('Rb3')
ax.title.set_text('Radius grid (um)')
plt.show()

#%%
idx = []

idx.append(np.where(rb_arr[:,0]==rb_arr[0,0]))
idx.append(np.where(rb_arr[:,1]==rb_arr[0,0]))
idx.append(np.where(rb_arr[:,2]==rb_arr[0,0]))

idx.append(np.where(rb_arr[:,0]==rb_arr[-1,-1]))
idx.append(np.where(rb_arr[:,1]==rb_arr[-1,-1]))
idx.append(np.where(rb_arr[:,2]==rb_arr[-1,-1]))

idx = np.array(idx).flatten()

idx = np.unique(idx)


rb_arr = rb_arr[idx,:]

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(rb_arr[:,0]*1e6,rb_arr[:,1]*1e6,rb_arr[:,2]*1e6, alpha=0.8,s=2)
ax.set_xlabel('Rb1')
ax.set_ylabel('Rb2')
ax.set_zlabel('Rb3')
ax.title.set_text('Radius grid (um)')
plt.show()
print(rb_arr.shape)

#%%

robot.verbose = False
for i in range(1):
    print(i)
    t = time.time()
    print(f' Inputs: {rb_arr[i,:]}')
    D1 = robot.forward_kinematics(rb_arr[i,0],rb_arr[i,1],rb_arr[i,2])
    print(f'elapsed time {i}/{rb_arr.shape[0]}: {(time.time()-t)*1e3:.2f} ms')



#%%
robot.save_data()     
#%%
robot.Plot_radius_grid()

#%%
robot.load_data()

[r_top, r_side1, r_side2] = robot.plot_workspace(ret_rad=True, l_threshold=65)


#%% Animate side view
idx_s1, idx_s2 = robot.get_idx_side_views(l_threshold=65)
data_plot = robot.data[idx_s1]
t = time.time()
robot.create_animation(data_plot,name='Side_WS_zoom',zoom=True, mesh=True)
print(f'elapsed time: {time.time()-t}')

#%% Animate Top view
idx_s1 = robot.get_idx_top_view(l_threshold=75)
data_plot = robot.data[idx_s1]
t = time.time()
robot.create_animation(data_plot,name='Top_WS_zoom',zoom=True, mesh=True)
print(f'elapsed time: {time.time()-t}')

#%%

data_plot = np.array(robot.data)
t = time.time()
robot.create_animation(data_plot,zoom=True, mesh=True)
print(f'elapsed time: {time.time()-t}')


#%%
robot.create_animation(data_plot,zoom=False)

#%%
[r_top, r_side1, r_side2] = robot.plot_workspace(ret_rad=True)

#%% Create path to test  Inverse Kinematics
L_sq = 20
h_plane = 0.5*(np.max(robot.End_effector[:,2])+np.min(robot.End_effector[:,2])) #970
h_plane = np.round(h_plane/10,0)*10
N = 20
x_offset = -380
N2 = int(N/2)
N4 = int(N/4)
st = np.linspace(-L_sq,L_sq,N,endpoint=True)
# --- Letter F
x = np.concatenate((np.ones(N)*st[0],st[:int(N/2)]))+L_sq*3/4
y = np.concatenate((st,np.ones(int(N/2))*st[-1]))
z = np.ones_like(x)*h_plane
r_inv = np.vstack((x,y,z)).T
##%%
x2 = np.concatenate((st[int(N/2)-1::-1],np.ones(int(N/2))*st[0],st[:int(N/4)]))+L_sq*3/4
y2 = np.concatenate((np.ones(int(N/2))*st[-1],st[:int(N/2)-1:-1], np.ones(int(N/4))*st[int(N/2)]))
z2 = np.ones_like(x)*h_plane

x_f = np.concatenate((x,x2))+x_offset
y_f = np.concatenate((y,y2))
z_f = np.ones_like(x_f)*h_plane
r_inv = np.vstack((x_f,y_f,z_f)).T

np.savetxt('F_letter_XYZ_input.txt', r_inv )

anim_xy(x_f,y_f)

#%% -- Letter E
x = np.concatenate((st[:int(N/2)][::-1],np.ones(N2)*st[0], st[:N4]))+L_sq*3/4
y = np.concatenate((np.ones(N2)*st[-1],st[:N2-1:-1], np.ones(N4)*st[N2]))
z = np.ones_like(x)*h_plane
r_inv = np.vstack((x,y,z)).T
#plt.scatter(r_inv[:,0],r_inv[:,1])

x2 = np.concatenate((st[:N4][::-1], np.ones(N2)*st[0], st[:N2]))+L_sq*3/4
y2 = np.concatenate((np.ones(N4)*st[N2], st[:N2][::-1],np.ones(N2)*st[0] ))
z2 = np.ones_like(x)*h_plane
x_e = np.concatenate((x,x2))
y_e = np.concatenate((y,y2))
z_e = np.ones_like(x)*h_plane
r_inv = np.vstack((x,y,z)).T
anim_xy(x_e,y_e)
#%% Letter M
x_m = np.concatenate((np.ones(N)*st[0], st[:N4],st[N4:2*N4],np.ones(N)*st[2*N4-1]))+L_sq/2
y_m = np.concatenate((st, st[-N4:][::-1],st[-N4:],st[::-1]))
z_m = np.ones_like(x)*h_plane
r_inv = np.vstack((x,y,z)).T
anim_xy(x_m,y_m)

#%% Letter T ---
x_t = np.concatenate((st[:N2],st[-N2-N4:-N2][::-1] , np.ones(N)*st[N4-1]))+L_sq/2
y_t = np.concatenate((np.ones(N2)*st[-1],np.ones(N4)*st[-1] , st[::-1]))
z_t = np.ones_like(x)*h_plane
r_inv = np.vstack((x,y,z)).T
anim_xy(x_t,y_t)
#%% Letter O ---
x_o = np.concatenate((st[:N2],np.ones(N)*st[N2],st[:N2][::-1],np.ones(N)*st[0]))+L_sq/2
y_o = np.concatenate((np.ones(N2)*st[0],st,np.ones(N2)*st[-1],st[::-1]))
z_o = np.ones_like(x)*h_plane
r_inv = np.vstack((x,y,z)).T
anim_xy(x_o,y_o)
#%%
x_tot = [x_f,x_e,x_m,x_t,x_0]
y_tot = [y_f,y_e,y_m,y_t,y_0]

#%%


# To save the animation using Pillow as a gif
# writer = animation.PillowWriter(fps=15,
#                                 metadata=dict(artist='Me'),
#                                 bitrate=1800)
# ani.save('scatter.gif', writer=writer)

#%% Solve paths
plt.scatter(r_inv[:,0],r_inv[:,1])
#plt.scatter(x,y)

#%% Solve Inverse kinematics of the points
s = []
for X_inv in r_inv[:]:
    t = time.time()
    S_inv = robot.Inverse_kinematics(X_inv)
    print(f'elapsed time: {time.time()-t}')
    s.append(S_inv)   
    print('Solved')

#%%

data_plot = np.array(s)

np.savetxt('F_letter_rads_05.txt', data_plot[:,-3:])
np.savetxt('F_letter_rads_SOLUTIONS5.txt', data_plot)

#data_plot = np.loadtxt('F_letter_rads_SOLUTIONS5.txt')

data_h = data_plot.copy()
s, i = robot.aproximate_HvsR()
def linear_interpolation(x, slope=s, intercept=i ):
    return slope * x + intercept

h1 = linear_interpolation(data_plot[:,-3]*1e6)*1e-6
h2 = linear_interpolation(data_plot[:,-2]*1e6)*1e-6
h3 = linear_interpolation(data_plot[:,-1]*1e6)*1e-6

data_h[:,-3] = h1
data_h[:,-2] = h2
data_h[:,-1] = h3
np.savetxt('F_letter_heights_05.txt', data_h[:,-3:])



#%%



#%%
plt.plot(data_plot[:,-3:])
plt.hlines(r_min*1e-6,0,55)
plt.hlines(r_max*1e-6,0,55)

#%%
t = time.time()
#robot.create_animation(data_plot,zoom=True, mesh=True)
robot.create_animation(data_plot,name='IK-Letter-F-Zoom', zoom=True, mesh=True)
print(f'elapsed time: {time.time()-t}')




#%% Graph r vs h


#%%%
from scipy.spatial import ConvexHull

points = robot.End_effector
hull = ConvexHull(points)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
#ax.scatter(robot.End_effector[:,0],robot.End_effector[:,1],robot.End_effector[:,2]) #points = self.End_effector[:,:2]
for simplex in hull.simplices:
    ax.plot3D(points[simplex, 0], points[simplex, 1], points[simplex, 2], 's-')
plt.show()



#%% ----- For the initial robot------
# -----------------------------------
# -----------------------------------


P_end = np.array([0,0,535])
robot = manip.Manip(P_end=P_end)
## Plot figures
[r_top, r_side1, r_side2] = robot.plot_workspace(ret_rad=True)

#%%
Size_md = 15*15*20
WS_md = 7.01
ration_md = WS_md/Size_md
Size_tp = 1.1*3*3
WS_tp = hull.volume*1e-9
ration_tp = WS_tp/Size_tp

#%%
figsize=(2,1.5)
fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111)
ax.plot(r_top, label=['r1','r2','r3'])
ax.legend()
f, (ax1, ax2) = plt.subplots(1, 2,figsize=(figsize[0]*2,figsize[1]), sharey=True)
ax1.plot(r_side1, label=['r1','r2','r3'])
ax1.legend()
ax2.plot(r_side2, label=['r1','r2','r3'])
ax2.legend()
#%% Solve points for forward kinematics
rb_l = np.arange(min(robot.data[:,6]), max(robot.data[:,6])+1e-7, 5e-6) 
rb_l = np.around(rb_l,6)
# Create mesh grid
rb_arr = np.array(np.meshgrid(rb_l,rb_l,rb_l)).T.reshape(-1,3)

rb_arr = np.array(np.meshgrid(rb_l,min(rb_l),rb_l)).T.reshape(-1,3)

print(rb_arr.shape)
# plot grid

pts = rb_arr

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(pts[:,0],pts[:,1],pts[:,2], alpha=0.8,s=2)
plt.show()


#%%

for i in range(rb_arr.shape[0]):
    
    t = time.time()
    D1 = robot.forward_kinematics(rb_arr[i,0],rb_arr[i,1],rb_arr[i,2])
    print(f'elapsed time {i}/{rb_arr.shape[0]}: {(time.time()-t)*1e3:.2f} ms')
    
#%%
robot.save_data()    

#%%
robot.plot_3D_workspace()



#%% plot grid

pts = robot.data[:,6:].copy()


fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(pts[:,0],pts[:,1],pts[:,2], alpha=0.1,s=2)
plt.show()

#%% Path inverse kinematics (Square)
L_sq = 40
N = 5
st = np.linspace(-L_sq,L_sq,N,endpoint=True)
x = np.concatenate((st,np.ones(N)*st[-1],st[::-1],np.ones(N)*st[0]))
y = np.concatenate((np.ones(N)*st[0],st,np.ones(N)*st[-1],st[::-1]))
z = np.ones_like(x)*1025
r_inv = np.vstack((x,y,z)).T
plt.scatter(r_inv[:,0],r_inv[:,1])
#plt.scatter(x,y)

#%% Solve Inverse kinematics of the points
s = []
for X_inv in r_inv[:]:
    t = time.time()
    S_inv = robot.Inverse_kinematics(X_inv)
    print(f'elapsed time: {time.time()-t}')
    s.append(S_inv)   
    print('Solved')
    
#%%    

data_plot = np.array(s)
t = time.time()
robot.create_animation(data_plot,zoom=True, mesh=False)
print(f'elapsed time: {time.time()-t}')

#%% Animate side view
idx_s1, idx_s2 = robot.get_idx_side_views()
data_plot = robot.data[idx_s1]
t = time.time()
robot.create_animation(data_plot,name='Side_WS_zoom',zoom=True, mesh=True)
print(f'elapsed time: {time.time()-t}')

#%% Animate Top view
idx_s1 = robot.get_idx_top_view()
data_plot = robot.data[idx_s1]
t = time.time()
robot.create_animation(data_plot,name='Top_WS_zoom',zoom=True, mesh=True)
print(f'elapsed time: {time.time()-t}')