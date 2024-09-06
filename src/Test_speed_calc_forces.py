#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 10:55:04 2024

@author: hector.ortiz
"""
from forces_eq import *

t = time.time()
x = 0
y = 0 
z = 388
ang_x = 0*np.pi/180
ang_y = 0*np.pi/180
ang_z = 0*np.pi/180
rb = 235e-6
F_l = []
#%%
t_l = []
for i in range(1):
    t = time.time()
    A = calc_forces(x, y, z, ang_x, ang_y, ang_z, rb)
    F_l.append(A)
    t_e = (time.time()-t)*1e3
    t_l.append(t_e)
    print(f'Elapsed time: {t_e} (ms)')
    r_err =((np.array(F_l[-1])-np.array(F_l[0]))/np.array(F_l[-1])*100).tolist()
    formatted_r_err = [f'{x:.2f}' for x in r_err]
    print(f'Rel. error %: {formatted_r_err}')
    
#%%
print(f'Mean Elapsed  time {np.mean(np.array(t_l))} ms')