#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 16:12:39 2024

@author: hector.ortiz
"""
import numpy as np
h = 480
r = 930.7/2
V = np.pi/6*h*(3*r**2+h**2)
print(f'Vol: {V*1e-9:0.3f} uL')