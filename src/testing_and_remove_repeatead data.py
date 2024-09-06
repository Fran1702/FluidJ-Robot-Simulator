#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 12:20:01 2024

@author: hector.ortiz
"""

# test 3d interpolation
from scipy.interpolate import LinearNDInterpolator

values = robot.data[:,:6]
x = robot.data[:,6]
y = robot.data[:,7]
z = robot.data[:,8]
#points = (x, y, z)
values = robot.data[:,:6]
points = robot.data[:,6:]
interp = LinearNDInterpolator(points, values, rescale=True)
#%%
point = robot.data[10,6:] + 1e-5
print(robot.data[10,:])
interp = LinearNDInterpolator(points, values)
print(interp(point))


#%%
arr = np.random.random((4,4,4,4))
x1 = np.array([0, 1, 2, 3])
x2 = np.array([0, 10, 20, 30])
x3 = np.array([0, 10, 20, 30])
x4 = np.array([0, .1, .2, .30])
points = (x1, x2, x3, x4)

#%%
point = robot.data[10,6:]
x_tuples = tuple(point)
y_tuples = [tuple(row) for row in points]

index = y_tuples.index(x_tuples)

#%%
# Remove repetead data and just keep the last calculated
inv_data = robot.data[:, 6:].copy()[::-1]
unique_values, unique_indices = np.unique(inv_data, axis=0, return_index=True)
unique_indices = inv_data.shape[0]-unique_indices-1
unique_values_last = robot.data[unique_indices]

robot.data = unique_values_last
# unique_values_last will contain the unique values with their last occurrences