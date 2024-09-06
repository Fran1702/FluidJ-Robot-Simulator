#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 10:37:54 2024

@author: hector.ortiz
"""
import forces_eq
import manip
import os
import numpy as np
#import pyqtgraph as pg
import time
import matplotlib.pyplot as plt

#%%
P_end = np.array([0,0,535])
robot = manip.Manip(P_end=P_end)
# Load Data
robot.load_data()
# Plot figures
robot.plot_workspace()
#%%
N = 10
X = robot.data[N][:6].copy()
R = robot.data[N][6:].copy()
R[0] = R[0]*1.2
robot.verbose = False
t = time.time()
S1 = robot.eqsystem_forward(X, *R )
print(f'elapsed time: {time.time()-t}')
F = S1[:3]
T = S1[3:]

#%% Testing Jacobian forward

S2 = robot.Jacob_forward(X, *R )

#%%
#Check forward kinematics:
    
robot.verbose = False

#robot.forward_kinematics(R[0],R[1],R[2], Q0=X.copy())
t = time.time()
D1 = robot.forward_kinematics( robot.data[10,6],robot.data[10,7],robot.data[10,8]*1.1)
print(f'elapsed time: {time.time()-t}')


#%% Test inverse

#Q is {rb{3},u(3)}
robot.verbose = False
U = robot.data[N][3:6].copy()
X_inv = np.concatenate((R, U))
Y_inv = robot.data[N][:3]
Y_inv = robot.End_effector[N]
t = time.time()
S2 = robot.eqsystem_inverse(X_inv, *Y_inv )
print(f'elapsed time: {time.time()-t}')
print(f'Eq system forward-inverse error: {S1-S2}')

#%% Test inverse kinematics
R_N = forces_eq.R_from_vec(D1[3:6].T)
P_N = D1[:3].copy()
X_inv = R_N@robot.P_end + P_N
#X_inv = robot.End_effector[0]
t = time.time()
S_inv = robot.Inverse_kinematics(X_inv)
print(f'elapsed time: {time.time()-t}')
print(f'forward-inverse error: {(D1-S_inv)}')

#%% Test using square
L_sq = 40
N = 5
st = np.linspace(-L_sq,L_sq,N,endpoint=True)
x = np.concatenate((st,np.ones(N)*st[-1],st[::-1],np.ones(N)*st[0]))
y = np.concatenate((np.ones(N)*st[0],st,np.ones(N)*st[-1],st[::-1]))
z = np.ones_like(x)*1025
r_inv = np.vstack((x,y,z)).T
plt.scatter(r_inv[:,0],r_inv[:,1])
#plt.scatter(x,y)

#%%
s = []
for X_inv in r_inv[:]:
    t = time.time()
    S_inv = robot.Inverse_kinematics(X_inv)
    print(f'elapsed time: {time.time()-t}')
    s.append(S_inv)   
    print('Solved')
    
#%%
robot.OUTPUT_FLAG = True
for i in range(len(s)) :
    R = np.array(s)[i,6:]
    X = np.array(s)[i,:6]
    print(robot.eqsystem_forward(X, *R))
    robot.save_mesh(i)
    
#%%    

data_plot = robot.data[:2,:]
#data_plot = np.vstack([data_plot,data_plot[0]])
#data_plot
t = time.time()
robot.create_animation(data_plot,zoom=True)
print(f'elapsed time: {time.time()-t}')



#%% Testing Jacobian inverse
#robot.OUTPUT_FLAG = False
S2 = robot.Jacob_inverse(X_inv, *Y_inv )
#%%
s2 = []
#for X_inv in r_inv[:]:
robot.verbose = False
S_inv = robot.Inverse_kinematics(r_inv[0], verbose=True)
s2.append(S_inv)   
print('Solved')
    

