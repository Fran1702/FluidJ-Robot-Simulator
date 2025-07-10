#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 09:16:26 2024

@author: hector.ortiz
"""
import numpy as np
from scipy.linalg import expm,logm
import subprocess
import os
import pandas as pd

import os
import sys

# Add the directory containing this file to the Python path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

OFF_SET_XY  = 500 # I add this offset to help the algorithm
U_OFFSET = 1
Z_OFFSET = 0
r_tripod = 850e-6
OUTPUT_FLAG = True
printing = False
ZMAX = 700
ZMIN = 200


def skew(x):
    return np.array([[0, -x[2], x[1]],
                     [x[2], 0, -x[0]],
                     [-x[1], x[0], 0]])


#%%

def theta_u_fromR(R):
    theta = np.arccos((np.trace(R)-1)/2)
    if theta!= 0:
        u = 1/(2*np.sin(theta))*np.array([R[2,1]-R[1,2], R[0,2]-R[2,0], R[1,0]-R[0,1]])
    else:
        u = np.array([0,0,0])
    return theta, u

#%%
def R_from_vec(u):
    # Calc rotation matrix using exponential map 
    s = np.shape(u)
    if len(s)>1:
        M = np.zeros((s[1],3,3))
        M[:,1,0] = u[2,:]
        M[:,0,1] = -u[2,:]
        M[:,2,0] = -u[1,:]
        M[:,0,2] = u[1,:]
        M[:,2,1] = u[0,:]
        M[:,1,2] = -u[0,:]
    else:
        M = np.zeros((3,3))
        M[1,0] = u[2]
        M[0,1] = -u[2]
        M[2,0] = -u[1]
        M[0,2] = u[1]
        M[2,1] = u[0]
        M[1,2] = -u[0]
    #print(M)

    
    #return expm(np.cross(np.eye(3),u))
    return expm(M)



def print_sol(sol):
    u_r = sol[0:3].copy()
    print('Ur: ', u_r)
    R = R_from_vec(u_r)
    phi, theta, psi = rot_2_angles(R)
    ang_x = phi
    ang_y = theta
    ang_z = psi
    txt = rf'$\phi_x$: {ang_x:.2f}$^\circ$, $\theta_y$: {ang_y:.2f}$^\circ$, $\psi_z$: {ang_z:.2f}$^\circ$'
    print(txt)
    P = sol[3:6].copy()
    P[0] = P[0]#-OFF_SET_XY
    P[1] = P[1]#-OFF_SET_XY

    txt = rf'$P_x$: {P[0]:.2f}$\mu m$, $P_y$: {P[1]:.2f}$\mu m$, $P_z$: {P[2]:.2f}$\mu m$'
    H1 = sol[6:9].copy()
    H1[1] = H1[1]#-OFF_SET_XY # Remove the offset
    H1[0] = H1[0]#-OFF_SET_XY # Remove the offset
    
    print('H1: ', H1)
    H2 = sol[9:12].copy()
    H2[1] = H2[1]#-OFF_SET_XY # Remove the offset
    H2[0] = H2[0]#-OFF_SET_XY # Remove the offset
    print('H2: ', H2)
    
    H3 = sol[12:15].copy()
    H3[1] = H3[1]#-OFF_SET_XY # Remove the offset
    H3[0] = H3[0]#-OFF_SET_XY # Remove the offset
    print('H3: ', H3)


def calc_R(phi, theta, psi):
    r11 = np.cos(psi)*np.cos(theta)
    r12 = np.cos(psi)*np.sin(theta)*np.sin(phi) - np.sin(psi)*np.cos(phi)
    r13 = np.cos(psi)*np.sin(theta)*np.cos(phi) + np.sin(psi)*np.sin(phi)
    
    r21 = np.sin(psi)*np.cos(theta)
    r22 = np.sin(psi)*np.sin(theta)*np.sin(phi) + np.cos(psi)*np.cos(phi)
    r23 = np.sin(psi)*np.sin(theta)*np.cos(phi) - np.cos(psi)*np.sin(phi)
    r31 = -np.sin(theta)
    r32 = np.cos(theta)*np.sin(phi)
    r33 = np.cos(theta)*np.cos(phi)
    
    R = np.array([[r11, r12, r13],
                 [r21, r22, r23],
                 [r31, r32, r33]])
    return R

def rot_2_angles(R, deg=True):
    
    theta = -np.arcsin(R[2,0])*180/np.pi
    
    phi = np.arctan2(R[2,1],R[2,2])*180/np.pi
    #if np.abs(R[2,0]) != 1:
    psi = np.arctan2(R[1,0]/np.cos(np.pi/180*theta),R[0,0]/np.cos(np.pi/180*theta))*180/np.pi
    #else:
    if deg:
        return phi, theta,psi 
    else:
        return phi*np.pi/180, theta*np.pi/180 , psi*np.pi/180 
    
def change_coord(x,y,z,ang_x, ang_y, ang_z, rb):
    with open('fluid_joint.fe', 'r') as file:
      content = file.readlines()
    
    x0 = x
    z0 = z
    y0 = y
    a_x = ang_x
    a_y = ang_y
    a_z = ang_z
    content[5] = f'PARAMETER RBOT = {rb:0.10e} \n'
    content[8] = f'PARAMETER ZMAX = {z0:0.10e} \n'
    content[9] = f'PARAMETER ang_y = {a_y:0.10e} \n'
    content[10] = f'PARAMETER ang_x = {a_x:0.10e} \n'
    content[11] = f'PARAMETER ang_z = {a_z:0.10e} \n'
    content[12] = f'parameter x0 = {x0:0.10e} // y-coord. \n'
    content[13] = f'parameter y0 = {y0:0.10e} // y-coord. \n'

    # and write everything back
    with open('fluid_joint.fe', 'w') as file:
        file.writelines(content)
        
def read_force():
    # Read the force file and return the last value
    
    fname = 'Forces.out'
    data = pd.read_csv(fname, sep=" ", header=None,skipinitialspace = True) #np.loadtxt(fname)
    
    fx = data.iloc[-1, 9]
    fy = data.iloc[-1, 10]
    fz = data.iloc[-1, 11]
    tx = data.iloc[-1, 12]
    ty = data.iloc[-1, 13]
    tz = data.iloc[-1, 14]
    
    '''    
    data = data[data.columns[:-1]]
    data = data.to_numpy()
    fy = data[-1,-4]
    fz = data[-1,-3]
    tx = data[-1,-2]
    ty = data[-1,-1]
    '''    
    
    
    return fx, fy, fz, tx, ty, tz

def check_calculated(y,z,ang_x, ang_y,rb):
    # Check if already exist the value, if exist return the forces, if not -1,-1
    
    
    fname = 'Forces.out'
    
    
    data = np.loadtxt(fname)
    rb_data = data[:,2]
    y_data = data[:,4]
    z_data = data[:,3]
    ang_x_data = data[:,5]
    ang_y_data = data[:,6]
    #print('y,z,rb, ang_x, ang_y',y,z,rb,ang_x, ang_y)
    
    '''
    y_data[y_data < 1e-10] = 0
    y = 0 if np.abs(y)< 1e-10 else y
    rb_str = "{:.7e}".format(rb)
    rb1, rb2 = rb_str.split('e')
    rb = float(rb1[:-1]+'e'+rb2)
    
    y_str = "{:.7e}".format(y)
    y1, y2 = y_str.split('e')
    y = float(y1[:-1]+'e'+y2)
    z_str = "{:.7e}".format(z)
    z1, z2 = z_str.split('e')
    z = float(z1[:-1]+'e'+z2)
    '''
    
    '''
    if sum((y_data==y)&(z_data==z)&(rb_data==rb)&(ang_x_data==ang_x)&(ang_y_data==ang_y))>0:
        #print('FOUND')
        #print('y, z, rb:', y, z, rb)
        fy = data[(y_data==y)&(z_data==z)&(rb_data==rb)&(ang_x_data==ang_x)&(ang_y_data==ang_y),-4][0]
        fz = data[(y_data==y)&(z_data==z)&(rb_data==rb)&(ang_x_data==ang_x)&(ang_y_data==ang_y),-3][0]
        tx = data[(y_data==y)&(z_data==z)&(rb_data==rb)&(ang_x_data==ang_x)&(ang_y_data==ang_y),-2][0]
        ty = data[(y_data==y)&(z_data==z)&(rb_data==rb)&(ang_x_data==ang_x)&(ang_y_data==ang_y),-1][0]
        #print('fy, fz, tx, ty: ', fy,fz, tx,ty)
        return fy,fz, tx, ty
    '''
    
    mask = (y_data == y) & (z_data == z) & (rb_data == rb) & (ang_x_data == ang_x) & (ang_y_data == ang_y)
    if mask.any():
        index = mask.argmax()
        fy, fz, tx, ty, tz = data[index, -5], data[index, -4], data[index, -3], data[index, -2], data[index, -1]
        return fy, fz, tx, ty, tz

    else:
        return -1,-1,-1,-1, -1

def calc_forces(x, y, z, ang_x, ang_y, ang_z, rb, outputfiles=False, j=None):
    # x,y,z in um
    # angles in rads
    #y = np.round(y,3)
    #z = np.round(z,3)
    
    #a_x = np.round(ang_x,6).astype(float)
    #a_x = a_x if np.abs(a_x)>1e-3 else 0.0
    
    #a_y = np.round(ang_y,6).astype(float)
    #a_y = a_y if np.abs(a_y)>1e-3 else 0.0
    
    #rb = np.round(rb*1e6,3)
    
#    if z>ZMAX:
#        return 0e10,-1e10, 0e10, 0e10, 0e10,0.0
#    if z<ZMIN:
#        return 0e10,1e10, 0e10, 0e10, 0e10, 0.0
    
    #change_coord(x*1e-6, y*1e-6,z*1e-6, ang_x, ang_y,ang_z, rb)
    
    #fy, fz,tx,ty,tz = check_calculated(y,z, a_x, a_y,rb)
    
    if outputfiles:
        # First check if the file exist:
        flag_before = os.path.isfile('vertex.txt')
        process = subprocess.Popen('evolver fluid_joint.fe',
                             stdin = subprocess.PIPE,
                             stdout = subprocess.PIPE, 
                             stderr = subprocess.PIPE,
                             text = True,
                             shell = True
                             )
        
        txt = f'RBOT := {rb:0.10e} \n recalc \n'
        txt = txt+f'new_z0 := {z*1e-6:0.10e} \n change_z0 \n'
        txt = txt+ f'new_ang_y := {ang_y:0.10e} \n change_ang_y \n'
        txt = txt+ f'new_ang_x := {ang_x:0.10e} \n change_ang_x \n'
        txt = txt+ f'new_ang_z := {ang_z:0.10e} \n change_ang_z \n'
        txt = txt+ f'new_x0 := {x*1e-6:0.10e} \n change_x0 \n'
        txt = txt+ f'new_y0 := {y*1e-6:0.10e} \n change_y0 \n'
        if j is not None:
            txt = txt+ f'vertex_fname := "vertex_{j}.txt" \n '
            txt = txt+ f'facet_fname := "facet_{j}.txt"  \n'
        inputdata = txt+"recalc \n run_forces \n save \n"
        #print(inputdata)
        #inputdata="save \n"
        stdoutdata,stderrdata=process.communicate(input=inputdata)
        #print(stdoutdata)
        lines = stdoutdata.split('\n')

        values = {
            'xforce': None,
            'yforce': None,
            'zforce': None,
            'xtorque': None,
            'ytorque': None,
            'ztorque': None
        }

        # Extract values by key
        for line in lines:
            for key in values:
                if line.startswith(f"{key}:"):
                    values[key] = float(line.split(':')[1].strip())

        # Make sure all values were found (optional but good practice)
        missing = [k for k, v in values.items() if v is None]
        if missing:
            raise ValueError(f"Missing values for: {', '.join(missing)}")

        # Now build the list in the desired order
        forces_l = [
            values['xforce'],
            values['yforce'],
            values['zforce'],
            values['xtorque'],
            values['ytorque'],
            values['ztorque'],
        ]
        #print(forces_l)
        return forces_l
        
    process = subprocess.Popen('evolver fluid_joint.fe',
                         stdin = subprocess.PIPE,
                         stdout = subprocess.PIPE, 
                         stderr = subprocess.PIPE,
                         text = True,
                         shell = True
                         )
    txt = f'RBOT := {rb:0.10e} \n recalc \n'
    txt = txt+f'new_z0 := {z*1e-6:0.10e} \n change_z0 \n'
    txt = txt+ f'new_ang_y := {ang_y:0.10e} \n change_ang_y \n'
    txt = txt+ f'new_ang_x := {ang_x:0.10e} \n change_ang_x \n'
    txt = txt+ f'new_ang_z := {ang_z:0.10e} \n change_ang_z \n'
    txt = txt+ f'new_x0 := {x*1e-6:0.10e} \n change_x0 \n'
    txt = txt+ f'new_y0 := {y*1e-6:0.10e} \n change_y0 \n'
    inputdata = txt+"recalc \n run_forces \n"
    
    #print(inputdata)
    stdoutdata,stderrdata=process.communicate(input=inputdata)
    
    l = stdoutdata.split('\n')
    #print([(l[i-7].split(' ')[-1]) for i in range(6)])
    try:
        forces_l = [float(l[i-7].split(' ')[-1]) for i in range(6)]
    except:
        print('Error computing forces:')
        print('INPUT DATA')
        print(inputdata)
        print('OUTPUT DATA')
        print(stdoutdata)
        
    for i in range(3):
        forces_l[i] = forces_l[i] if np.abs(forces_l[i])>1e-8 else 0
    for i in range(3):
        forces_l[i+3] = forces_l[i+3] if np.abs(forces_l[i+3])>1e-12 else 0
        
    #print(forces_l)
    #fx, fy, fz, tx, ty, tz = read_force()
#    input()
    return forces_l

#%%
#%%

def eqsystem(Q,*args):
    rbase1,rbase2,rbase3 = args
    if printing:
        print(Q)
        input()
    #input()
    Q = Q.copy()
    # Q is {u{3},P{3}}
    # rb is the base radius of each droplet
    #u_vec = np.ones(3)
    #u_vec[0:2] = Q[0:2]
    u_vec = Q[3:]
    u_vec = u_vec - U_OFFSET# scale the values to improve solution
    R = R_from_vec(u_vec)
    a_x, a_y, a_z = rot_2_angles(R, deg=False)
    P = Q[0:3]
    P[0] = P[0]-OFF_SET_XY
    P[1] = P[1]-OFF_SET_XY
    P[2] = P[2]+Z_OFFSET # Remove the offset
    
    w1 = np.array([r_tripod,0,0])*1e6
    w2 = np.array([-0.5, np.sqrt(3)/2, 0])*r_tripod*1e6   # um
    w3 = np.array([-0.5, -np.sqrt(3)/2, 0])*r_tripod*1e6
    b1 = w1.copy()
    b2 = w2.copy()
    b3 = w3.copy()
    
    # Kinematics equation [um]
    H1 = R@b1 + P - w1
    H2 = R@b2 + P - w2
    H3 = R@b3 + P - w3
    if printing:
        print('H1: ',H1)
        print('H2: ',H2)
        print('H3: ',H3)
    
    '''
    K1 = (H1+w1-R@b1-P)
    K2 = (H2+w2-R@b2-P)
    K3 = (H3+w3-R@b3-P)
    
    if printing:
        print('K1: ',K1)
        print('K2: ',K2)
        print('K3: ',K3)
        
    H1 = Q[6:9]
    H1[1] = H1[1]-OFF_SET_XY # Remove the offset
    H1[0] = H1[0]-OFF_SET_XY # Remove the offset
    H1[2] = H1[2]+Z_OFFSET # Remove the offset
    H2 = Q[9:12]
    H2[1] = H2[1]-OFF_SET_XY # Remove the offset
    H2[0] = H2[0]-OFF_SET_XY # Remove the offset
    H2[2] = H2[2]+Z_OFFSET # Remove the offset
    H3 = Q[12:15]
    H3[1] = H3[1]-OFF_SET_XY # Remove the offset
    H3[0] = H3[0]-OFF_SET_XY # Remove the offset
    H3[2] = H3[2]+Z_OFFSET # Remove the offset
    '''

    
    # First calculate the forces:
    # Transform the xy displacement to one displacement as the axis are rotated
    # To calculate the torques is needed the normal vectors and its angles
    # u_0 = R*u_t+P0 ; u_0 vector in fixed reference frame
    # u_0 = R_se*u_se ; u_t vector in tripod frame, u_se in se frame
    # u_se = R*R_se.T*u_t
    # u_t = R.T*R_se*u_se
    #
    # (only displacement in y axis of the droplet is allowed)
    f = []
    t = []
    n_tripod = np.cross((R@b1-R@b2),(R@b3-R@b2)) # Normal vector to the plane of the tripod (x_t, y_x)
    
    
    for (H, rbase) in zip([H1,H2,H3],[rbase1,rbase2,rbase3]):
        x = H[0]
        y = H[1]
        z = H[2]
        d = np.sqrt(x**2+y**2)
        d_vec = np.array([x,y,0])
        if d < 1e-10:
            d=0
        # Angle between the vector displacement and the y axis of fixed coordinate system
        if d!= 0:
            #angle_se = np.arccos(np.dot(np.array([0,1,0]),d_vec)/d)
            v1 = d_vec/d
            v2 = np.array([0,1,0])
            R_ang = np.array([[v1[0]*v2[0]+v1[1]*v2[1], v2[0]*v1[1]-v1[0]*v2[1],0],
                           [v1[0]*v2[1]-v2[0]*v1[1], v1[0]*v2[0]+v1[1]*v2[1],0],
                           [0,0,1]])
        else:
            angle_se = 0
            R_ang = R_from_vec(angle_se*np.array([0,0,1]))
            if printing:
                print('D is zero')
        # Rotatation matrix  to rotate in z axis 

        R_se = R@(R_ang.T)
        
        phi, theta,psi = rot_2_angles(R_se)
        ang_x = phi 
        ang_y = theta
        ang_z = psi
        if printing:
            txt = rf'$\phi_x$: {ang_x:.2f}$^\circ$, $\theta_y$: {ang_y:.2f}$^\circ$, $\psi_z$: {ang_z:.2f}$^\circ$'
            print(txt)
        ang_x = phi*np.pi/180
        ang_y = theta*np.pi/180
        ang_z = psi*np.pi/180
        #print(' X,Y, Z, rb, ang_x, ang_y VAL: ', x,y,z, rb, ang_x, ang_y)
        #input()
#        if z != 1050.0 or z!= 1100.0:
            #input()

        #fy, fz, tx_se, ty_se, tz_se = calc_forces(d,z, ang_x, ang_y,rbase, outputfiles=OUTPUT_FLAG)
        fx, fy, fz, tx, ty, tz = calc_forces(x,y,z, a_x, a_y ,a_z, rbase, outputfiles=OUTPUT_FLAG)
        # projection to original xy axis
        '''
        ft = np.sqrt(fy**2+fz**2)
        if d !=0:
            #fx = np.abs(x/d)*fy
            #fy = np.abs(y/d)*fy
            #tx = np.abs(x/d)*ty_se - np.abs(y/d)*tx_se
            #ty = np.abs(y/d)*ty_se + np.abs(x/d)*tx_se
            #tz = tz_se
            f_fixed = R_se@np.array([0,fy,fz])
            t_fixed = R_se@np.array([tx_se,ty_se,tz_se])
            #print('f fixed', f_fixed)
            #print('tau fixed', t_fixed)
            #print(f'Tx: {tx}, Ty: {ty}')
            #print(f'Tx_se: {tx_se}, Ty_se: {ty_se}')
            
        else:
            f_fixed = np.array([0,0,fz])
            t_fixed = np.array([0,0,0])
        #f.append(f_fixed)
        #t.append(t_fixed)
        '''
        f.append(np.array([fx,fy,fz]))
        
        t.append(np.array([tx,ty,tz]))
        #if printing:
         #   input('i')
        
    # Force equation:
    F = f[0]+f[1]+f[2] # [N]
    #print(f[0])
    #print(f[1])
    #print(f[2])
    F = F*1e6 # To uN
    if printing:
        print('F: ',F)
    
    # Torque equations f in N b in um:
    T = np.cross(R@b1, f[0]) + np.cross(R@b2, f[1]) + np.cross(R@b3, f[2])
    T = T+t[0]*1e6+t[1]*1e6+t[2]*1e6
    if printing:
        print('Torq vecs 1: ',np.cross(R@b1+P, f[0]))
        print('Torq vecs 2: ',np.cross(R@b2+P, f[1]) )
        print('Torq vecs 3: ', np.cross(R@b3+P, f[2]))
        print('t1,t2,t3', t[0], t[1], t[2])
        print('T: ',T)
        
    T = T # [uNm]
    

    
    # Rotation equation
    #R1 = R@(R.T)-np.identity(3)
#    if printing:
#        print('R:', R1.flatten())


    #print(f'Norm: {np.linalg.norm(np.concatenate((K1,K2,K3,T,F)).flatten()):.2e}')
    #print(f'Norm: {np.linalg.norm(np.concatenate((T,F)).flatten()):.2e}')
    #return np.linalg.norm(np.concatenate((K1,K2,K3,T,F)).flatten())
    #return np.concatenate((K1,K2,K3,T,F)).flatten()
    #return np.concatenate((F,T,K1,K2,K3)).flatten()
    #return np.concatenate((F,T[:2])).flatten()
    return np.concatenate((F,T[:3])).flatten()

def eqsystem_inverse(Q,*args):
    # Q is {rb{3},u(3)}
    # args is P
    #rbase1,rbase2,rbase3 = args
    rbase1,rbase2,rbase3 = Q[3:]
    P = args
    if printing:
        print(Q)
        input()
    #input()
    Q = Q.copy()
    # Q is {u{3},P{3}}
    # rb is the base radius of each droplet
    #u_vec = np.ones(3)
    #u_vec[0:2] = Q[0:2]
    u_vec = Q[3:]
    u_vec = u_vec - U_OFFSET# scale the values to improve solution
    R = R_from_vec(u_vec)
    a_x, a_y, a_z = rot_2_angles(R, deg=False)
    P = Q[0:3]
    P[0] = P[0]-OFF_SET_XY
    P[1] = P[1]-OFF_SET_XY
    P[2] = P[2]+Z_OFFSET # Remove the offset
    
    w1 = np.array([r_tripod,0,0])*1e6
    w2 = np.array([-0.5, np.sqrt(3)/2, 0])*r_tripod*1e6   # um
    w3 = np.array([-0.5, -np.sqrt(3)/2, 0])*r_tripod*1e6
    b1 = w1.copy()
    b2 = w2.copy()
    b3 = w3.copy()
    
    # Kinematics equation [um]
    H1 = R@b1 + P - w1
    H2 = R@b2 + P - w2
    H3 = R@b3 + P - w3
    if printing:
        print('H1: ',H1)
        print('H2: ',H2)
        print('H3: ',H3)
    
    '''
    K1 = (H1+w1-R@b1-P)
    K2 = (H2+w2-R@b2-P)
    K3 = (H3+w3-R@b3-P)
    
    if printing:
        print('K1: ',K1)
        print('K2: ',K2)
        print('K3: ',K3)
        
    H1 = Q[6:9]
    H1[1] = H1[1]-OFF_SET_XY # Remove the offset
    H1[0] = H1[0]-OFF_SET_XY # Remove the offset
    H1[2] = H1[2]+Z_OFFSET # Remove the offset
    H2 = Q[9:12]
    H2[1] = H2[1]-OFF_SET_XY # Remove the offset
    H2[0] = H2[0]-OFF_SET_XY # Remove the offset
    H2[2] = H2[2]+Z_OFFSET # Remove the offset
    H3 = Q[12:15]
    H3[1] = H3[1]-OFF_SET_XY # Remove the offset
    H3[0] = H3[0]-OFF_SET_XY # Remove the offset
    H3[2] = H3[2]+Z_OFFSET # Remove the offset
    '''

    
    # First calculate the forces:
    # Transform the xy displacement to one displacement as the axis are rotated
    # To calculate the torques is needed the normal vectors and its angles
    # u_0 = R*u_t+P0 ; u_0 vector in fixed reference frame
    # u_0 = R_se*u_se ; u_t vector in tripod frame, u_se in se frame
    # u_se = R*R_se.T*u_t
    # u_t = R.T*R_se*u_se
    #
    # (only displacement in y axis of the droplet is allowed)
    f = []
    t = []
    n_tripod = np.cross((R@b1-R@b2),(R@b3-R@b2)) # Normal vector to the plane of the tripod (x_t, y_x)
    
    
    for (H, rbase) in zip([H1,H2,H3],[rbase1,rbase2,rbase3]):
        x = H[0]
        y = H[1]
        z = H[2]
        d = np.sqrt(x**2+y**2)
        d_vec = np.array([x,y,0])
        if d < 1e-10:
            d=0
        # Angle between the vector displacement and the y axis of fixed coordinate system
        if d!= 0:
            #angle_se = np.arccos(np.dot(np.array([0,1,0]),d_vec)/d)
            v1 = d_vec/d
            v2 = np.array([0,1,0])
            R_ang = np.array([[v1[0]*v2[0]+v1[1]*v2[1], v2[0]*v1[1]-v1[0]*v2[1],0],
                           [v1[0]*v2[1]-v2[0]*v1[1], v1[0]*v2[0]+v1[1]*v2[1],0],
                           [0,0,1]])
        else:
            angle_se = 0
            R_ang = R_from_vec(angle_se*np.array([0,0,1]))
            if printing:
                print('D is zero')
        # Rotatation matrix  to rotate in z axis 

        R_se = R@(R_ang.T)
        
        phi, theta,psi = rot_2_angles(R_se)
        ang_x = phi 
        ang_y = theta
        ang_z = psi
        if printing:
            txt = rf'$\phi_x$: {ang_x:.2f}$^\circ$, $\theta_y$: {ang_y:.2f}$^\circ$, $\psi_z$: {ang_z:.2f}$^\circ$'
            print(txt)
        ang_x = phi*np.pi/180
        ang_y = theta*np.pi/180
        ang_z = psi*np.pi/180
        #print(' X,Y, Z, rb, ang_x, ang_y VAL: ', x,y,z, rb, ang_x, ang_y)
        #input()
#        if z != 1050.0 or z!= 1100.0:
            #input()

        #fy, fz, tx_se, ty_se, tz_se = calc_forces(d,z, ang_x, ang_y,rbase, outputfiles=OUTPUT_FLAG)
        fx, fy, fz, tx, ty, tz = calc_forces(x,y,z, a_x, a_y ,a_z, rbase, outputfiles=OUTPUT_FLAG)
        # projection to original xy axis
        '''
        ft = np.sqrt(fy**2+fz**2)
        if d !=0:
            #fx = np.abs(x/d)*fy
            #fy = np.abs(y/d)*fy
            #tx = np.abs(x/d)*ty_se - np.abs(y/d)*tx_se
            #ty = np.abs(y/d)*ty_se + np.abs(x/d)*tx_se
            #tz = tz_se
            f_fixed = R_se@np.array([0,fy,fz])
            t_fixed = R_se@np.array([tx_se,ty_se,tz_se])
            #print('f fixed', f_fixed)
            #print('tau fixed', t_fixed)
            #print(f'Tx: {tx}, Ty: {ty}')
            #print(f'Tx_se: {tx_se}, Ty_se: {ty_se}')
            
        else:
            f_fixed = np.array([0,0,fz])
            t_fixed = np.array([0,0,0])
        #f.append(f_fixed)
        #t.append(t_fixed)
        '''
        f.append(np.array([fx,fy,fz]))
        
        t.append(np.array([tx,ty,tz]))
        #if printing:
         #   input('i')
        
    # Force equation:
    F = f[0]+f[1]+f[2] # [N]
    #print(f[0])
    #print(f[1])
    #print(f[2])
    F = F*1e6 # To uN
    if printing:
        print('F: ',F)
    
    # Torque equations f in N b in um:
    T = np.cross(R@b1, f[0]) + np.cross(R@b2, f[1]) + np.cross(R@b3, f[2])
    T = T+t[0]*1e6+t[1]*1e6+t[2]*1e6
    if printing:
        print('Torq vecs 1: ',np.cross(R@b1+P, f[0]))
        print('Torq vecs 2: ',np.cross(R@b2+P, f[1]) )
        print('Torq vecs 3: ', np.cross(R@b3+P, f[2]))
        print('t1,t2,t3', t[0], t[1], t[2])
        print('T: ',T)
        
    T = T # [uNm]
    

    
    # Rotation equation
    #R1 = R@(R.T)-np.identity(3)
#    if printing:
#        print('R:', R1.flatten())


    #print(f'Norm: {np.linalg.norm(np.concatenate((K1,K2,K3,T,F)).flatten()):.2e}')
    #print(f'Norm: {np.linalg.norm(np.concatenate((T,F)).flatten()):.2e}')
    #return np.linalg.norm(np.concatenate((K1,K2,K3,T,F)).flatten())
    #return np.concatenate((K1,K2,K3,T,F)).flatten()
    #return np.concatenate((F,T,K1,K2,K3)).flatten()
    #return np.concatenate((F,T[:2])).flatten()
    return np.concatenate((F,T[:3])).flatten()
