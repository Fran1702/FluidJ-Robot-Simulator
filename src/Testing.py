#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 12:34:18 2024

@author: hector.ortiz
"""
from concave_hull import concave_hull, concave_hull_indexes

#%%
points = robot.End_effector[:,::2]
idxes = concave_hull_indexes(
    points,
    length_threshold=25,
)

#%%

rb_l = np.arange(robot.R_DROPLET_MIN, robot.R_DROPLET_MAX+1e-6, 10e-6) 
rb1 = np.concatenate((rb_l, np.ones(len(rb_l)*2)*rb_l[0],
                      rb_l,np.ones(len(rb_l)*2)*rb_l[-1]))
rb2 = np.concatenate((np.ones(len(rb_l)*1)*rb_l[0], rb_l[::-1],
                      np.ones(len(rb_l)*2)*rb_l[0], rb_l,
                      np.ones(len(rb_l)*1)*rb_l[-1]))

rb_arr = np.vstack((rb1,rb2,rb3)).T # Just if used generated before
x = np.around(rb_arr,6)
y = np.around(res[:,6:],6)
idx= []
for i in range(x.shape[0]):
    ind = np.all(y[:] == x[i],axis=1).nonzero()[0]
    if len(ind)>0:
        idx.append(ind[0])
        
idx_side = np.array(idx).flatten()  

#%%
plt.scatter(robot.End_effector[idxes,0],robot.End_effector[idxes,2],s=4)
plt.scatter(robot.End_effector[:,0],robot.End_effector[:,2], s=1)

#%%
plt.plot(robot.data[idxes,6])
plt.plot(robot.data[idxes,7])
plt.plot(robot.data[idxes,8])