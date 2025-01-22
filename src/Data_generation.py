#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 09:16:32 2024

@author: hector.ortiz
"""

import forces_eq
import manip
import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import alphashape 
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.colors import LightSource

from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

import pyvista as pv
import pandas as pd

#%%
P_end = np.array([-380,0,630]) # ANTENNA as end effector
#P_end = np.array([0,0,0])
#DROPLET_VOLUME = 0.12e-9#0.22e-9
#DROPLET_VOLUME = 0.22e-9
DROPLET_VOLUME = 0.102e-9#0.22e-9
#DROPLET_VOLUME = 0.3e-9#0.22e-9
theta_0 = 135*np.pi/180 # 160
theta_min = 60*np.pi/180 # 65
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
rb_l = np.linspace(robot.R_DROPLET_MIN*1e-6, robot.R_DROPLET_MAX*1e-6,19)
#rb_l = np.linspace(robot.R_DROPLET_MIN*1e-6, robot.R_DROPLET_MAX*1e-6,11)
rb_l = np.around(rb_l,6)
# Create mesh grid
rb_arr = np.array(np.meshgrid(rb_l,rb_l,rb_l)).T.reshape(-1,3)

#rb_arr = np.array(np.meshgrid(rb_l,min(rb_l),rb_l)).T.reshape(-1,3)

print(rb_arr.shape)

pts = rb_arr

#fig = plt.figure()
#ax = fig.add_subplot(111, projection="3d")
#ax.scatter(rb_arr[:,0]*1e6,rb_arr[:,1]*1e6,rb_arr[:,2]*1e6, alpha=0.8,s=2)
#ax.set_xlabel('Rb1')
#ax.set_ylabel('Rb2')
#ax.set_zlabel('Rb3')
#ax.title.set_text('Radius grid (um)')
#plt.show()

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

#fig = plt.figure()
#ax = fig.add_subplot(111, projection="3d")
#ax.scatter(rb_arr[:,0]*1e6,rb_arr[:,1]*1e6,rb_arr[:,2]*1e6, alpha=0.8,s=2)
#ax.set_xlabel('Rb1')
#ax.set_ylabel('Rb2')
#ax.set_zlabel('Rb3')
#ax.title.set_text('Radius grid (um)')
#plt.show()
print(rb_arr.shape)



#%%
robot.verbose = False
# Assuming rb_arr is a numpy array
shuffled_indices = np.random.permutation(rb_arr.shape[0])
j=0
for i in shuffled_indices:
#for i in range(rb_arr.shape[0]):
#for i in range(10):
    #i = 378
    print(f'Solving: {j}/{rb_arr.shape[0]}')
    t = time.time()
    print(f'Inputs: {rb_arr[i,:]}')
    D1 = robot.forward_kinematics(rb_arr[i,0],rb_arr[i,1],rb_arr[i,2])
    print(f'elapsed time {j}/{rb_arr.shape[0]}: {(time.time()-t)*1e3:.2f} ms')
    j = j+1
    robot.save_data()         


