#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 11:13:36 2024

@author: hector.ortiz
"""

#import forces_eq
import manip
import numpy as np
import time
#import os
#import matplotlib.pyplot as plt
#import matplotlib.animation as animation

#import alphashape 
#from mpl_toolkits.mplot3d import Axes3D
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#from mpl_toolkits.mplot3d.art3d import Line3DCollection
#from matplotlib.colors import LightSource

#from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
#from matplotlib.gridspec import GridSpec
#import pyvista as pv
#import pandas as pd
#from Closed_loop import *
#from Plotting_CL import *
import random


#%% ------------- Robot definition and data load  -----------------
#P_end = np.array([0,0,630])
P_end = np.array([-380,0,630]) # ANTENNA as end effector
DROPLET_VOLUME = 0.102e-9#0.22e-9
theta_0 = 135*np.pi/180
theta_min = 60*np.pi/180
r_min = (3*DROPLET_VOLUME*np.sin(theta_0)**2/(np.pi*(2+np.cos(theta_0))*(1-np.cos(theta_0))**2))**(1/3)*1e6
r_max = (3*DROPLET_VOLUME*np.sin(theta_min)**2/(np.pi*(2+np.cos(theta_min))*(1-np.cos(theta_min))**2))**(1/3)*1e6
r_min = np.round(r_min*1.0,0)
r_max = np.round(r_max*1.052,0)
print(r_min)
print(r_max)
robot = manip.Manip(P_end=P_end, V0=DROPLET_VOLUME,R_DROPLET_MIN=r_min, R_DROPLET_MAX=r_max,                    
                    ZMAX=600, ZMIN=200, R_top=150)
# Load Data
robot.load_data()
#%%
rb_arr = robot.data[:,-3:]
N_data = rb_arr.shape[0]
indices = list(range(N_data))
random.shuffle(indices)  # Shuffle the indices to ensure all are processed exactly once
k = -1
for i in indices:  # Iterate over the shuffled indices
#for i in range(10):
    #i = 378
    k=k+1
    print(f'Processing data:{k}/{N_data}, index number: {i}')    
    t = time.time()
    print(f' Inputs: {rb_arr[i,:]}')
    D1 = robot.forward_kinematics(rb_arr[i,0],rb_arr[i,1],rb_arr[i,2])
    print(f'elapsed time {k}/{rb_arr.shape[0]}: {(time.time()-t)*1e3:.2f} ms')
    robot.save_data()         