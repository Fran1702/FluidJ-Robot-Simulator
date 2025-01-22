#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 12:58:53 2024

@author: fran
"""

import sys

# Add the directory containing this file to the Python path


from forces_eq import *
import numpy as np
import time
from mpl_toolkits.mplot3d import Axes3D, art3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from stl import mesh
from mpl_toolkits import mplot3d
from scipy.optimize import fsolve, root_scalar
import matplotlib.ticker as ticker
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from concave_hull import concave_hull, concave_hull_indexes
np.set_printoptions(suppress=True, formatter={'float': '{:0.5e}'.format})
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
from multiprocessing.pool import ThreadPool
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from scipy.interpolate import NearestNDInterpolator
from math import log10, floor
os.chdir(os.path.dirname(os.path.abspath(__file__)))
#%%
import pyvista as pv
from pyvistaqt import BackgroundPlotter
from pv_plotter_extra import *
#%%
import scienceplots
plt.style.use(['science','ieee'])
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
#plt.rcParams.update({'figure.dpi':'100'})

mpl.rcParams['legend.frameon'] = 'True'
mpl.rcParams['legend.facecolor'] = 'w'

# For 3d plot
import WS_3D
import mpl_toolkits.mplot3d as a3

def get_mean_change(data, axis):
    # return mean and max deviation
    l1 = np.max(data[:,axis])
    l2 = np.min(data[:,axis])
    return (l1+l2)/2, l1-l2

def matprint(mat, fmt="g"):
    col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
        print("")
  
# Try it!      

class Manip:
    
    def __init__(self, OFF_SET_XY  = 0,  U_OFFSET = 1, Z_OFFSET = 0, 
                 r_tripod = 850e-6, OUTPUT_FLAG = False, verbose = False,
                 ZMAX = 700, ZMIN = 200, R_DROPLET_MIN=488, R_DROPLET_MAX=588,
                 V0=0.3e-9, R_top = 150, P_end = [0,0,0],
                 delete_data =False):
        # Asign variables and flags
        self.OFF_SET_XY = OFF_SET_XY
        self.OFF_SET_XY_SOL = 0
        self.U_OFFSET = U_OFFSET
        self.U_OFFSET_SOL = 0
        self.Z_OFFSET = Z_OFFSET
        self.Z_OFFSET_SOL = 0
        self.r_tripod = r_tripod
        self.OUTPUT_FLAG = OUTPUT_FLAG
        self.verbose = verbose
        self.R_DROPLET_MIN = R_DROPLET_MIN
        self.R_DROPLET_MAX = R_DROPLET_MAX
        self.V0 = V0
        self.R_top = R_top # um
        self.ZMAX = ZMAX
        self.ZMIN = ZMIN
        self.delete_data = delete_data
        self.w1 = np.array([self.r_tripod,0,0])*1e6
        self.w2 = np.array([-0.5, np.sqrt(3)/2, 0])*self.r_tripod*1e6   # um
        self.w3 = np.array([-0.5, -np.sqrt(3)/2, 0])*self.r_tripod*1e6
        self.b1 = self.w1.copy()
        self.b2 = self.w2.copy()
        self.b3 = self.w3.copy()
        self.tol_solved = 5e-3
        self.P_end = P_end # Coordinates of End efector in movile coordinates
        self.data = None
        self.Init_SE_file() # init SE file with the Volume value
        self.Init_LinearModel()
    
    def Init_SE_file(self):
        #sys.path.append(os.path.dirname(os.path.abspath(__file__)))
        
        #print(os.path.dirname(os.path.abspath(__file__)))
        # Add  an update of the PARAMETER ZMAX to avoid explotion of the 
        # surface It can be the radius of a half sphere of volume V...
        with open('fluid_joint.fe', 'r') as file:
            content = file.readlines()
        # Set volume in SE file
        a = ["PARAMETER V0" in s for s in content]
        idx = np.where(a)[0][0] # Index where V0 is located
        content[idx] = f'PARAMETER V0 = {self.V0:0.3e}  // mm3 \n'
        # Set Rtop in SE file
        a = ["PARAMETER RTOP" in s for s in content]
        idx = np.where(a)[0][0] # Index where V0 is located
        content[idx] = f'PARAMETER RTOP = {self.R_top*1e-6:0.3e}  // \n'
        # Set Z in SE file
        a = ["PARAMETER ZMAX" in s for s in content]
        idx = np.where(a)[0][0] # Index where V0 is located
        content[idx] = f'PARAMETER ZMAX = {self.ZMAX*1e-6:0.3e}  // \n'
        
        
        
        # and write everything back
        with open('fluid_joint.fe', 'w') as file:
            file.writelines(content)
            
            
        
        return
          
    def Plot_radius_grid(self):
        
        rb_arr = self.data[:,-3:]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(rb_arr[:,0]*1e6,rb_arr[:,1]*1e6,rb_arr[:,2]*1e6, alpha=0.8,s=2)
        ax.set_xlabel('Rb1')
        ax.set_ylabel('Rb2')
        ax.set_zlabel('Rb3')
        ax.title.set_text('Radius grid (um)')
        plt.show()

    def aproximate_HvsR(self):
        '''
         Linear approximation of the height vs radius of the droplet
        
        '''
        # Extract the last three columns
        last_three_cols = self.data[:, -3:]

        # Check where all values in the last three columns are equal
        equal_indices = np.all(last_three_cols == last_three_cols[:, [0]], axis=1)

        # Get the indices where the condition is met
        indices = np.where(equal_indices)[0]
        r = self.data[indices,-3]*1e6
        z = self.data[indices, 2]
        # Perform linear regression to find the slope and intercept
        slope, intercept = np.polyfit(r, z, 1)

        # Generate the y-values of the best-fit line
        y_fit = slope * r + intercept
        # Plot the best-fit line
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(r,z, '.',label='Data Points')
        ax.plot(r, y_fit, '--' ,label='Best-fit Line')
        
        # Annotate slope and intercept on the graph
        annotation_text = f'Equation: y = {slope:.2f}x + {intercept:.2f}'
        plt.text(0.05, 0.15, annotation_text, transform=plt.gca().transAxes, 
                 fontsize=8, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.6))
        
        plt.legend()
        plt.show()
        
        return slope, intercept
        
    def Init_LinearModel(self):
        # Implement this funtion to run an initilization of the linear model
        #name = f'data_F_T_{self'
        #np.savetxt('../data/Solutions_WS.txt',self.data)
        return 
    
    def load_data(self):
        ''' Load the Data with shape: 
            P_1x3, U_1x3, r1, r2,r3 shape: 9,1'''
            
        #os.chdir("../data")
        try:
            self.data = np.loadtxt(f'../data/Solutions_WS_{self.V0*1e9}uL.txt')
            # Delete data out of the range
            if self.delete_data:
                for i in range(3):
                    idx = np.where(self.data[:,-(i+1)] >  self.R_DROPLET_MAX*1.001*1e-6)    
                
                    self.data = np.delete(self.data, idx, 0) # Remove data
                    idx = np.where(self.data[:,-(i+1)] < self.R_DROPLET_MIN*0.999*1e-6)    
                    self.data = np.delete(self.data, idx, 0) # Remove data
                

            R_N = R_from_vec(self.data[:,3:6].T)
            P_N = self.data[:,:3].copy()
            self.End_effector = R_N@self.P_end + P_N
        except:
            
            print('No previus data of the model')
        return
    
    def save_mesh(self,i):      
        # Directory paths
        new_dir = "../data/mesh/"
        
        # Lists of base file names to be renamed
        files_to_rename = [
            ('facet', 3),
            ('vertex', 3)
        ]
        
        # Loop over files and perform batch renaming
        for file_prefix, count in files_to_rename:
            for j in range(count):
                old_file = f'{file_prefix}_{j+1}.txt'
                new_file = os.path.join(new_dir, f'{file_prefix}_{j+1}_{i}.txt')
                if os.path.exists(new_file):
                    # Delete the existing file
                    os.remove(new_file)
                os.rename(old_file, new_file)
        return   
    
    def get_idx_top_view(self,l_threshold=10):
        points = self.End_effector[:,:2]
        idx_top = concave_hull_indexes(
            points,
            length_threshold=l_threshold,
        )
        idx_top = np.append(np.array(idx_top),idx_top[0] )
        return idx_top
    
    def get_idx_side_views(self,l_threshold=10):
        points = self.End_effector[:,1:]
        idx_side = concave_hull_indexes(
            points,
            length_threshold=l_threshold,
        )
        idx_side = np.append(np.array(idx_side),idx_side[0] )
        points = self.End_effector[:,::2]
        idx_side2 = concave_hull_indexes(
            points,
            length_threshold=l_threshold,
        )
        idx_side2 = np.append(np.array(idx_side),idx_side[0] )
        return idx_side, idx_side2
    
    def plot_workspace(self, l_threshold=10,figsize=(2,1.5), ret_arg = False):
        '''
        Function that plot the workspace of the manipulator

        Parameters
        ----------
        l_threshold : TYPE, optional
            DESCRIPTION. The default is 10.
        figsize : TYPE, optional
            DESCRIPTION. The default is (2,1.5).
        ret_rad : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        list
            DESCRIPTION.

        '''
        # l_threshold, threshold for concave hull
        # Top view Plot
        points = self.End_effector[:,:2]
        idx_top = concave_hull_indexes(
            points,
            length_threshold=l_threshold,
        )
        idx_top = np.append(np.array(idx_top),idx_top[0] )
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        ax.plot(self.End_effector[idx_top,0], self.End_effector[idx_top,1],'-', alpha=1)
        ax.fill(self.End_effector[idx_top,0], self.End_effector[idx_top,1], color = 'gray', alpha = 0.3)
        ax.set_ylabel(r'Y ($\mu$m)')
        ax.set_xlabel(r'X ($\mu$m)')
        ax.set_title('Top view')
        plt.savefig("../data/figs/Top_WS.pdf", format="pdf", bbox_inches="tight")
        
        # Side view
        points = self.End_effector[:,1:]
        idx_side = concave_hull_indexes(
            points,
            length_threshold=l_threshold,
        )
        idx_side_YZ = np.append(np.array(idx_side),idx_side[0] )
        f, (ax1, ax2) = plt.subplots(1, 2,figsize=(figsize[0]*2,figsize[1]), sharey=True)
        ax2.plot(self.End_effector[idx_side_YZ,1],self.End_effector[idx_side_YZ,2],'-',zorder=20, alpha=1)
        ax2.fill(self.End_effector[idx_side_YZ,1],self.End_effector[idx_side_YZ,2], color = 'gray', alpha = 0.3)
        
        ax2.set_xlabel(r'Y ($\mu$m)')
        points = self.End_effector[:,::2]
        idx_side = concave_hull_indexes(
            points,
            length_threshold=l_threshold,
        )
        
        idx_side_XZ = np.append(np.array(idx_side),idx_side[0] )
        
        ax1.plot(self.End_effector[idx_side_XZ,0],self.End_effector[idx_side_XZ,2],'-',zorder=20, alpha=1)
        ax1.fill(self.End_effector[idx_side_XZ,0],self.End_effector[idx_side_XZ,2], color = 'gray', alpha = 0.3)

        ax1.set_ylabel(r'Z ($\mu$m)')
        ax1.set_xlabel(r'X ($\mu$m)')
        f.suptitle('Side view')
        plt.savefig("../data/figs/Side_WS.pdf", format="pdf", bbox_inches="tight")
        plt.show()
        if ret_arg==1:
            return [self.data[idx_top,-3:], self.data[idx_side_XZ,-3:],self.data[idx_side_YZ,-3:]]
        if ret_arg==2:
            return [idx_top, idx_side_XZ, idx_side_YZ]
        else:
            ret
        
    def eqsystem_thread(self,H, rbase,j=None):
        # rbase in m
        # x,y,z in um
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
                if self.verbose:
                    print('D is zero')
            # Rotatation matrix  to rotate in z axis 

            R_se = self.R@(R_ang.T)
            a_x, a_y, a_z = rot_2_angles(self.R, deg=False)
            if self.verbose:
                txt = rf'φx: {a_x*180/np.pi:.2f}$^\circ$, θy: {a_y*180/np.pi:.2f}$^\circ$, ψz: {a_z*180/np.pi:.2f}$^\circ$'
                print(txt)
            # Check extremes of Z 
            if z>self.ZMAX:
                z = self.ZMAX
            elif z<self.ZMIN:
                z = self.ZMIN
            # Check extremes of rbase
            # Here I just run some checks, to be improved with good value
            '''
            if rbase > 1.2*self.ZMAX*1e-6:
                print("rbase to high")
                rbase = 1.2*self.ZMAX*1e-6
                
            elif rbase < self.ZMIN*1e-6:
                print("rbase to low")
                rbase = self.ZMIN*1e-6
            '''
            fx, fy, fz, tx, ty, tz = calc_forces(x,y,z, a_x, a_y ,a_z, rbase, outputfiles=self.OUTPUT_FLAG, j=j)
            # projection to original xy axis
            f = np.array([fx,fy,fz])
            t = np.array([tx,ty,tz])
            return f, t
        
    def eqsystem_forward(self, Q, *args):
        # Receives the Q with the OFFSETS
        rbase1,rbase2,rbase3 = args
        if self.verbose:
            print(Q)
            
        print_time = False
        if print_time:
            t0 = time.time() 
        #print('Q', Q)    
          #  input()
        #input()
        Q = Q.copy()
        # Q is {P{3}, u{3}}
        # rb is the base radius of each droplet
        #u_vec = np.ones(3)
        #u_vec[0:2] = Q[0:2]
        u_vec = Q[3:]
        u_vec = u_vec - self.U_OFFSET_SOL  # remove scale the values to improve solution
        self.R = R_from_vec(u_vec)
        a_x, a_y, a_z = rot_2_angles(self.R, deg=False)
        P = Q[0:3]
        P[0] = P[0]-self.OFF_SET_XY_SOL
        P[1] = P[1]-self.OFF_SET_XY_SOL
        P[2] = P[2]-self.Z_OFFSET_SOL # Remove the offset
        

        
        # Kinematics equation [um]
        H1 = self.R@self.b1 + P - self.w1
        H2 = self.R@self.b2 + P - self.w2
        H3 = self.R@self.b3 + P - self.w3
        if self.verbose:
            print('H1: ',H1)
            print('H2: ',H2)
            print('H3: ',H3)
        
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
        n_tripod = np.cross((self.R@self.b1-self.R@self.b2),
                            (self.R@self.b3-self.R@self.b2)) # Normal vector to the plane of the tripod (x_t, y_x)
        
        # The following cycle can be improved using threads...
    #    for (H, rbase) in zip([H1,H2,H3],[rbase1,rbase2,rbase3]):
        with ThreadPool() as pool:
            args = [(H, rbase,j) for H, rbase,j in zip([H1, H2, H3], [rbase1, rbase2, rbase3],[1,2,3])]
            async_result = pool.starmap_async(self.eqsystem_thread, args)
            results = async_result.get()
        for result in results:
            #print(result)
            f.append(result[0])
            t.append(result[1])
        #if self.verbose:
         #   input('i')
        
        # Force equation:
        F = f[0]+f[1]+f[2] # [N]
        F = F*1e6 # To uN
        if self.verbose:
            print('F: ',F)
        
        # Torque equations f in N b in um:
        T = np.cross(self.R@self.b1, f[0]) + np.cross(self.R@self.b2, f[1]) + np.cross(self.R@self.b3, f[2])
        T = T+t[0]*1e6+t[1]*1e6+t[2]*1e6
        if self.verbose:
            print('Torq vecs 1: ',np.cross(self.R@self.b1+P, f[0]))
            print('Torq vecs 2: ',np.cross(self.R@self.b2+P, f[1]) )
            print('Torq vecs 3: ', np.cross(self.R@self.b3+P, f[2]))
            print('t1,t2,t3', t[0], t[1], t[2])
            print('T: ',T)
            
        T = T # [uNm]
        #print('F,T', np.concatenate((F,T[:3])).flatten())
        if print_time:
            print(f'Eq. system computed in {time.time()-t0:0.1f} (s)')
        return np.concatenate((F,T[:3])).flatten()
    
    def init_guess(self,rb1,rb2,rb3):
        debug = False
        h_tilde = []
        r_input = [rb1,rb2,rb3]
        if debug: print('R_input: ', r_input)
        for i in range(3):
            f_h = lambda x: -x**3-x*(3*r_input[i]**2 + 3*(self.R_top*1e-6)**2) + 6*self.V0/np.pi
            #print(root_scalar(f_h, x0=(self.ZMIN*1e-6+ self.ZMAX*1e-6)/2,bracket=[self.ZMIN*1e-6, self.ZMAX*1e-6]))
            
            h_approx = fsolve(f_h, (self.ZMIN*1e-6+ self.ZMAX*1e-6)/2)[0]*1e6
            h_tilde.append(h_approx*np.array([0,0,1]))
            
            #Z_approx = fsolve(f_h, (self.ZMIN*1e-6+ self.ZMAX*1e-6)/2)[0]*0.8*1e6
        if debug: print('h_tilde: ', h_tilde)
        v1 = self.b1 + h_tilde[0]
        v2 = self.b2 + h_tilde[1]
        v3 = self.b3 + h_tilde[2]
        P_tilde = 1/3*(v1+v2+v3)
        if debug:  print('P_Disp: ', P_tilde)
        z_tilde = np.cross(v3-v2,v1-v2)
        z_tilde = z_tilde/np.linalg.norm(z_tilde)
        if debug: print('z_tilde: ',z_tilde)
        v = np.cross(np.array([0,0,1]),z_tilde)
        if debug: print('v', v)
        s = np.linalg.norm(v)
        if debug: print('s', s)
        c = np.dot(np.array([0,0,1]), z_tilde)
        if debug: print('c',c)
        v_x = skew(v)
        if debug: print(v_x)
        R_tilde = np.eye(3) + v_x + v_x@v_x*(1/(1+c))
        if debug: print(R_tilde)
        t_tilde, u_tilde = theta_u_fromR(R_tilde)
        u_tilde = t_tilde*u_tilde
        if debug: print('t', t_tilde)
        if debug: print('u',u_tilde)
        
        P_tilde = R_tilde@P_tilde
        if debug:  print('P_tilde: ', P_tilde)
        Q0 = np.array([P_tilde[0], P_tilde[1],P_tilde[2], u_tilde[0]+self.U_OFFSET, u_tilde[1]+self.U_OFFSET, u_tilde[2]+self.U_OFFSET])
        
        return Q0
    
    
    
    def forward_kinematics(self,rb1,rb2,rb3, Q0=None):
        # Q is the position and orientation vector of a initial point[P,u]
        # Return the Position and orientation of the frame refrence mobile P and u
        self.OFF_SET_XY_SOL = self.OFF_SET_XY
        self.U_OFFSET_SOL = self.U_OFFSET
        self.Z_OFFSET_SOL = self.Z_OFFSET
        FLAG_SOLVE = True
        FLAG_SOL_INDATA = False
        FLAG_INTEPOLATED = False
        rb1 = round(rb1, -7*int(floor(log10(abs(rb1)))))
        rb2 = round(rb2, -7*int(floor(log10(abs(rb2)))))
        rb3 = round(rb3, -7*int(floor(log10(abs(rb3)))))
        rb_arr = np.array([rb1,rb2,rb3])
        if Q0 is not None:
            Q0[3:] = Q0[3:]+self.U_OFFSET_SOL
            Q0[0:2] = Q0[0:2] + self.OFF_SET_XY_SOL
            Q0[2] = Q0[2] +  self.Z_OFFSET_SOL
            
        else:
            # If res is defined, find the closest value and use it as a Q0
            # if not stimate one initial value (can fail)
            try: 
                
                if self.data is None:
                    raise
                    
                Q0_arr = self.data[:,:6].copy()
                Q0_arr[:,3:] = Q0_arr[:,3:]+self.U_OFFSET_SOL
                Q0_arr[:,0:2] = Q0_arr[:,0:2] + self.OFF_SET_XY
                Q0_arr[:,2] = Q0_arr[:,2] +  self.Z_OFFSET
                
                #Q0_in = self.data[:,6:]
                Q0_in = self.data[:,6:9]

                point = np.array([rb1,rb2,rb3])
                # Round radius
                #print('rb_arr',rb_arr)
                array1_reshaped = rb_arr[np.newaxis, :] 
                array2_reshaped = Q0_in[np.newaxis, :, :]  
    
                # Calculate the norms between each element of array1 and array2
                norms = np.linalg.norm(array1_reshaped - array2_reshaped, axis=2)
               # print('norms shape: ',norms.shape)
                # Find the index of the minimum norm for each element of array1
                closest_indices = np.argmin(norms, axis=1)  # Shape: (128,)
                norms = norms.T
               # print('closest_indices', closest_indices[0])
               # print('norm',norms[closest_indices[0]])
                #print('rb: ', array1_reshaped[closest_indices[0]])
                #print('data: ', array2_reshaped[closest_indices[0]])      
                Q0 = Q0_arr[closest_indices[0]].copy()
                # Check if the solution exist, return it forwardly and not continue
                #print('Q0: ',Q0)
                
                points = self.data[:,6:9].copy()
                #points = self.data[:,6:]
                
                # Check if the solution exist, return it forwardly and not continue
                try: 
                    
                    tolerance = 1e-6  # Example tolerance value
                    #print('norm',norms[closest_indices[0]])
                    if norms[closest_indices[0]]>tolerance:
                        #print(rb_arr)
                        #print(array2_reshaped[closest_indices[0]])
                        #print(norms[closest_indices[0]])
                        #print('FAR')
                        raise
                    index = closest_indices[0]
                    #index = y_tuples.index(x_tuples)
#                    x_tuples = tuple(point)
#                    print('x_tuples', x_tuples)
#                    y_tuples = [tuple(row) for row in points]
                    #print('y_tuples', y_tuples)
#                    index = y_tuples.index(x_tuples)
                    #print('ZXC')
                    Q0 = Q0_arr[index].copy()
                    #if index != 419:
                    #return self.data[index].copy() # Remove ths line to calculate it each time
                    FT = self.eqsystem_forward(Q0, *rb_arr)
                    if self.verbose:
                        print(f'Norm FT: {np.linalg.norm(FT)} ')
                    
                    if np.linalg.norm(FT) <= self.tol_solved:
                        print('Sol. in data')
                        FLAG_SOLVE = False
                    
                    if np.linalg.norm(FT) > self.tol_solved:
                        # Recalc
                        print('High norm')
                        #print(index)
                        self.data = np.delete(self.data, index, 0) # Remove bad data
                        raise # Exception to interpolet guess
                        
                    if self.verbose:
                        print('Sol. in data')
                        FLAG_SOL_INDATA = True
                        # Check if the forces are small
                        
                            
                    #return self.data[index].copy()
                except:
                # Linear interpolation
                    #FT = self.eqsystem_forward(Q0.copy(), *rb_arr)
                    
                    values =  self.data[:,:6].copy()
                    points = self.data[:,6:9].copy()
                    Q0 = self.init_guess(rb1,rb2,rb3)
                    #print('Q0 guess: ',Q0)
                    interp = LinearNDInterpolator(points, values, rescale=True)
                    #interp = RBFInterpolator(points, values, kernel='linear')
                    #print(interp(point))
                    Q0 = interp(point)[0].copy()
                    Q0[3:] = Q0[3:] + self.U_OFFSET
                    print('Guess interpolated')
                   # print('Q0 interp: ',Q0)
                    #print(Q0)
                    if any(np.isnan(Q0)):
                        print('NAN: Q0 aproximated: Using nearest and spherical')
                        interp = NearestNDInterpolator(points, values)
                        #print(interp(point))
                        Q0 = interp(point)[0].copy()
                        Q0[3:] = Q0[3:] + self.U_OFFSET        
                        Q0_sp = self.init_guess(rb1,rb2,rb3)   
                        Q0 = (Q0 + Q0_sp)/2
                        if self.verbose: print('Q0 guess: ',Q0)
                        if any(np.isnan(Q0)):
                            if self.verbose: print('NAN: Q0 aproximated: Using spherical cap ')
                        #print('Nearest')
                        # if is outside the convex hull, takes the mean between closest value
                        # And the using spherical cap approximation
                        #Q0 = Q0_arr[closest_indices[0]]
                            Q0 = self.init_guess(rb1,rb2,rb3)
                            #print('Q0 guess: ',Q0)
                    # check if the approximation is good enough
                    FT = self.eqsystem_forward(Q0, *rb_arr)                            
                    if np.linalg.norm(FT) <= self.tol_solved:
                        print('Interpolation is good')
                        FLAG_SOLVE = False
                        FLAG_INTEPOLATED = True
                        sol = Q0.copy()
                        sol[3:] = sol[3:] - self.U_OFFSET_SOL
                        sol[0:2] = sol[0:2] - self.OFF_SET_XY
                        output = np.concatenate((sol,np.array([rb1,rb2,rb3])))
                        #print(f'output: {output}')
                    else:
                        print(f'Interpolation not good enough: {np.linalg.norm(FT)}')
                 #   print(Q0)
            except:
                print('Q0 aproximated: Using nearest and s; 1')
                Q0_sp = self.init_guess(rb1,rb2,rb3)   
                #print(Q0_sp)
                #print('checkpoint 0') 
                #print(self.data)
                if self.data is not None:
                    
                    values =  self.data[:,:6].copy()
                    points = self.data[:,6:9]
                    interp = NearestNDInterpolator(points, values)
                    #print(interp(point))
                    Q0 = interp(point)[0].copy()
                    Q0[3:] = Q0[3:] + self.U_OFFSET        
                    Q0 = (Q0 + Q0_sp)/2
                else:
                    Q0 = Q0_sp.copy()
                
                #print('checkpoint 0')    
                
                if any(np.isnan(Q0)):
                    print('NAN: Q0 aproximated: Using spherical cap ')
                    #print('Q0 aproximated: Using spherical cap')
                    #print(np.mean([rb1,rb2,rb3]))
                    Q0 = self.init_guess(rb1,rb2,rb3) # With offset
                    Q0_clean = Q0.copy()
                    Q0_clean[3:] = Q0_clean[3:]-np.ones_like(Q0_clean[3:])*self.U_OFFSET
                    #print('Q0 guess: ',Q0_clean)
                    
                print('checkpoint 1')    
                FT = self.eqsystem_forward(Q0, *rb_arr)   
                
                if np.linalg.norm(FT) <= self.tol_solved:
                    print('Interpolation is good')
                    FLAG_SOLVE = False
                    FLAG_INTEPOLATED = True
                    sol = Q0.copy()
                    sol[3:] = sol[3:] - self.U_OFFSET_SOL
                    sol[0:2] = sol[0:2] - self.OFF_SET_XY
                    output = np.concatenate((sol,np.array([rb1,rb2,rb3])))
                    #print(f'output: {output}')
                    
                else:
                    print(f'Interpolation not good enough: {np.linalg.norm(FT)}')
                 
        if self.verbose:
            t0 = time.time() 
        
            
        if FLAG_SOLVE:
            ## Here use linear interpolation model to get a good first guess Q0
            # Before solve it, check if the interpolation is good:
                
            #FT = self.eqsystem_forward(Q0.copy(), *rb_arr)
            
                    
            res = fsolve(self.eqsystem_forward, Q0.copy(),args=tuple([rb1,rb2,rb3]), fprime=self.Jacob_forward,
                         epsfcn=1e-4, factor=100, full_output=True, xtol=1e-6)
        
            sol = res[0].copy()
            info = res
            if self.verbose:
                print(info)
                print(f'Norm: {np.linalg.norm(info[1]["fvec"]):0.2e}')
        # Remove offset of the solution
            sol[3:] = sol[3:] - self.U_OFFSET_SOL
            sol[0:2] = sol[0:2] - self.OFF_SET_XY
            if self.verbose:
                print(f'Solved in {time.time()-t0:0.1f} (s)')
            
        
        self.OFF_SET_XY_SOL = 0
        self.U_OFFSET_SOL = 0
        self.Z_OFFSET_SOL = 0
        
        
        # When is solved, add the data to the current data
        if FLAG_SOLVE or FLAG_INTEPOLATED:
            output = np.concatenate((sol,np.array([rb1,rb2,rb3])))
            if FLAG_SOLVE:
                print(f'Norm: {np.linalg.norm(info[1]["fvec"]):0.2e}')
            else:
                print(f'Norm: {np.linalg.norm(FT):0.2e}') 
        else:
            output = self.data[index].copy()
        #print(output.shape)
        if self.data is not None:
            # just save if was solved otherwise the sol exist in the database
            if FLAG_SOLVE or FLAG_INTEPOLATED:
                if FLAG_SOLVE:
                    if np.linalg.norm(info[1]["fvec"]) <  self.tol_solved: 
                        self.data = np.append(self.data,np.reshape(output,(1,9)),axis=0)
                        self.save_data()
                    else:
                        print('Norm to high not saved')
                elif FLAG_INTEPOLATED :
                    self.data = np.append(self.data,np.reshape(output,(1,9)),axis=0)
                
                
        else:
            print('HERE')
            self.data = np.reshape(output,(1,9))
        
        # Add the solution to the end effector list
        #self.End_effector = R_N@self.P_end + P_N   
        if not FLAG_SOL_INDATA or FLAG_INTEPOLATED:
            u_r = output[3:].copy()
            R = R_from_vec(u_r)
            P = output[0:3].copy()
            End_effector = R@self.P_end + P
            #End_effector = End_effector[:,np.newaxis]
            #print(End_effector.shape)
            if hasattr(self, 'End_effector'):
                self.End_effector = np.vstack([self.End_effector,End_effector.reshape(1, 3)])
            else:
                self.End_effector = np.vstack([End_effector.reshape(1, 3)])
                
        return output
 
    
    def eqsystem_inverse(self, Q, *args):
        
        # Q is {rb{3},u(3)}
        # args is r_end in fixed coordnates
        Px,Py,Pz = args
        rbase1,rbase2,rbase3 = Q[:3]
        r_end = np.array([Px,Py,Pz])
        
        if self.verbose:
            print(Q)
            print(rbase1,rbase2,rbase3)
          #  input()
        #input()
        #Q = Q.copy()
        # Q is {u{3},P{3}}
        # rb is the base radius of each droplet
        #u_vec = np.ones(3)
        #u_vec[0:2] = Q[0:2]
        u_vec = Q[3:]
        u_vec = u_vec - self.U_OFFSET_SOL  # scale the values to improve solution
        self.R = R_from_vec(u_vec)
        a_x, a_y, a_z = rot_2_angles(self.R, deg=False)
        # transform P_end to the translation of the frame reference P
        P = r_end - self.R@self.P_end

    
        # Kinematics equation [um]
        H1 = self.R@self.b1 + P - self.w1
        H2 = self.R@self.b2 + P - self.w2
        H3 = self.R@self.b3 + P - self.w3
        if self.verbose:
            print('H1: ',H1)
            print('H2: ',H2)
            print('H3: ',H3)
        
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
        n_tripod = np.cross((self.R@self.b1-self.R@self.b2),
                            (self.R@self.b3-self.R@self.b2)) # Normal vector to the plane of the tripod (x_t, y_x)
        
        # The following cycle can be improved using threads...
        with ThreadPool() as pool:
            args = [(H, rbase) for H, rbase in zip([H1, H2, H3], [rbase1, rbase2, rbase3])]
            async_result = pool.starmap_async(self.eqsystem_thread, args)
            results = async_result.get()
        for result in results:
            #print(result)

            f.append(result[0])
            t.append(result[1])
            
        # Force equation:
        F = f[0]+f[1]+f[2] # [N]
        F = F*1e6 # To uN
        if self.verbose:
            print('F: ',F)
        
        # Torque equations f in N b in um:
        T = np.cross(self.R@self.b1, f[0]) + np.cross(self.R@self.b2, f[1]) +\
            np.cross(self.R@self.b3, f[2])
        T = T+t[0]*1e6+t[1]*1e6+t[2]*1e6
        if self.verbose:
            print('Torq vecs 1: ',np.cross(self.R@self.b1+P, f[0]))
            print('Torq vecs 2: ',np.cross(self.R@self.b2+P, f[1]) )
            print('Torq vecs 3: ', np.cross(self.R@self.b3+P, f[2]))
            print('t1,t2,t3', t[0], t[1], t[2])
            print('T: ',T)
            
        T = T # [uNm]
        return np.concatenate((F,T[:3])).flatten()

    def Jacob_inverse(self,Q, *args):
        # Computes the jacobian for the inverse kinematics
        flag_debug = True
        if flag_debug: print("Jacob:")
        if flag_debug: print(Q)
        d = 2e-5 # radii in m
        d_ang = 0.1*np.pi/180.0 # angles
        JM = np.zeros((6,6))
        rx, ry, rz = args
        # create the thread pool
        diag = [d,d,d,d_ang,d_ang,d_ang]
        d_m = np.diag(diag)
        
        Q_l = np.array([Q.copy() for i in range(6)]) + d_m
        Q_l = np.append(Q_l,[Q.copy()],axis=0)
        with ThreadPool() as pool:
            args = [(Q_l[i],*args) for i in range(len(Q_l))]
            async_result = pool.starmap_async(self.eqsystem_inverse, args)
            results = async_result.get()
        # iterate return values and report
        dF = []
        for result in results:
            dF.append(result)    
        dF_l = []
        for i in range(6):
            dF_l.append(dF[i]-dF[-1])
        
        for i  in range(3):
            JM[:,i] = dF_l[i]/d
        for i  in range(3):
            JM[:,i+3] = dF_l[i+3]/d_ang
            
        if flag_debug: print('JAcobian: ')
        if flag_debug:  matprint(JM)
        
        return JM
    
    def Jacob_forward(self,Q, *args):
        # Computes the jacobian for the inverse kinematics
        # Receives the Q with OFFSETS
        flag_debug = False
        if flag_debug: print("Jacob:")
        if flag_debug: print(Q)
        d = 0.1 # um
        d_ang = 0.1*np.pi/180.0
        diag = [d,d,d,d_ang,d_ang,d_ang]
        d_m = np.diag(diag)
        JM = np.zeros((6,6))
        rx, ry, rz = args
        # create the thread pool
        print_time = False
        if print_time:
            t0 = time.time() 
        
        Q_l = np.array([Q.copy() for i in range(6)]) + d_m
        Q_l = np.append(Q_l,[Q.copy()],axis=0)
        with ThreadPool() as pool:
            args_pool = [(Q_l[i],*args) for i in range(len(Q_l))]
            async_result = pool.starmap_async(self.eqsystem_forward, args_pool)
            results = async_result.get()
        # iterate return values and report
        dF = []
        for result in results:
          #  print('results: ', result)
            dF.append(result)    
        
        for i in range(6):
            dF.append(dF[i]-dF[-1])
        
        for i  in range(3):
            JM[:,i] = dF[i]/d
        for i  in range(3):
            JM[:,i+3] = dF[i+3]/d_ang
        
        if print_time:
            print(f'Jacobian computed in {time.time()-t0:0.1f} (s)')
        if flag_debug: print('JAcobian: ')
        if flag_debug:  matprint(JM)
        return JM
    
    def Inverse_kinematics(self, r_end, X0=None, verbose=False):
        """Calculates the inverse kinematics

        Parameters
        ----------
        r_end : np.array
            Coordinates point of the end effector to calculates its inverse
        
        Returns
        -------
        np.array
        an array containing [P, u, rb1, rb2, rb3 ]
            """
        
        # r_end is the desire position of the END EFECTOR in fixed coordinates!
        
        # Q is {rb{3},u(3)}
        # args is P
        # Return the P
        self.U_OFFSET_SOL = self.U_OFFSET
        
        if X0 is not None:
            X0[3:] = X0[3:]+self.U_OFFSET_SOL
        
        else:
            # If res is defined, find the closest value and use it as a Q0
            # if not stimate one initial value (can fail)
            #try:   
            Q0_arr = self.data[:,np.array([6,7,8,3,4,5])].copy()
            Q0_arr[:,3:] = Q0_arr[:,3:]+self.U_OFFSET_SOL
            
            Q0_in = self.End_effector[:,:3]
            # Reshape arrays to enable broadcasting
            r_end_arr = np.array([r_end[0],r_end[1],r_end[2]])
            
            array1_reshaped = r_end_arr[np.newaxis, :] 
            array2_reshaped = Q0_in[np.newaxis, :, :]  

            # Calculate the norms between each element of array1 and array2
            norms = np.linalg.norm(array1_reshaped - array2_reshaped, axis=2)
            # Find the index of the minimum norm for each element of array1
            closest_indices = np.argmin(norms, axis=1)  # Shape: (128,)
            X0 = Q0_arr[closest_indices[0]]
            # Check if the solution exist, return it forwardly and not continue
            if norms[0,closest_indices[0]]==0:
                #np.concatenate((sol,np.array([rb1,rb2,rb3])))
                if self.verbose or verbose:
                    print('Sol. in data')
                return self.data[[closest_indices[0]]].copy()
            #except:
                #print('Q0 aproximated')
                #Q0 = np.array([self.OFF_SET_XY, self.OFF_SET_XY,(self.ZMAX+self.ZMIN)/2, 1,1,1])
        #return
    
        if self.verbose or verbose:
            t0 = time.time() 
            
            
        res = fsolve(self.eqsystem_inverse, X0.copy(), fprime=self.Jacob_inverse,
                     args=tuple([r_end[0],r_end[1],r_end[2]]),
                      epsfcn=1e-5, factor=100, full_output=True, xtol=1e-8)
        
        sol = res[0].copy()
        info = res
        if self.verbose or verbose:
            print(info)
            print(f'Norm: {np.linalg.norm(info[1]["fvec"]):0.2e}')
        # Remove offset of the solution
        sol[3:] = sol[3:]-1
        if self.verbose or verbose:
            print(f'Solved in {time.time()-t0:0.1f} (s)')
        
        
        self.U_OFFSET_SOL = 0
        u = sol[3:]
        P = r_end-R_from_vec(u)@self.P_end
        
        return np.concatenate((P,sol[3:],sol[:3]))

    def init(self):
        pass
    
    def animate(self, i, data_plot):
        if self.verbose:
            print('Frame: ', i)
        resf = data_plot[i]
        
        self.ax.cla()
        self.ax2.cla()
        self.ax3.cla()
        #u_r = resf[0:3].copy()
        u_r = resf[3:].copy()
        if self.verbose:
            print('Ur: ', u_r)
        R = R_from_vec(u_r)
        phi, theta, psi = rot_2_angles(R)
        ang_x = phi
        ang_y = theta
        ang_z = psi
        
        txt = rf'$\phi_x$: {ang_x:.2f}$^\circ$, $\theta_y$: {ang_y:.2f}$^\circ$, $\psi_z$: {ang_z:.2f}$^\circ$'
        if self.verbose:
            print(txt)
        P = resf[0:3].copy()
        txt = txt+'\n'+rf'$P_x$: {P[0]:.2f} $\mu m$, $P_y$: {P[1]:.2f} $\mu m$, $P_z$: {P[2]:.2f} $\mu m$'
        self.ax.text(-1.5e3,-1.5e3, -1.9e3, txt, 'x', fontsize=7,
                bbox=dict(facecolor='white', alpha=1, edgecolor='black', linewidth=0.2, boxstyle='square'))
        # Kinematics equation [um]
        H1 = R@self.b1 + P - self.w1
        H2 = R@self.b2 + P - self.w2
        H3 = R@self.b3 + P - self.w3
        if self.verbose:
            print(f'P = {P}')
            print(f'H1 = {H1}')
        #P_end = np.array([0,0,535]) # 100 um in z 
        if self.verbose:
            print(f'Pend = {P_end}')
        End_effector = R@self.P_end + P
        self.P_limits.append(End_effector)
        if self.verbose:
            print(f'End Eff = {End_effector}')
        theta = np.linspace(0, 2 * np.pi, 201)
        j = 0
        for (v2,v1) in zip([H1, H2, H3], [self.w1,self.w2,self.w3]):
            j = j+1
            l = np.linalg.norm(v1-v2)
           
        ### Plot droplets
           
            file_path = f'../data/mesh/vertex_{j}_{i}.txt'
            data = pd.read_csv(file_path, sep=' ', header=None).to_numpy()
            
            idx_v = data[:,0].astype(int)
            idx_n = np.arange(len(idx_v))
            #verts = data[:,1:]*1e6
            verts = data[:,:].copy()
            verts[:,1:] = verts[:,1:]*1e6
            
            
            verts[:,1] = verts[:,1]+v1[0]
            verts[:,2] = verts[:,2]+v1[1]
            verts_dict = dict(zip(verts[:, 0], verts[:, 1:]))
            
            data = np.loadtxt(f'../data/mesh/facet_{j}_{i}.txt')
            
            facets = data.copy()
            faces_orig = data[:,1:4].astype(int)
            cols = data[:,-1].astype(int)
            faces = faces_orig.copy()
    
            val = ['blue','red','black', 'white']

            faces = faces[cols[:].argsort()]
            cols = cols[cols[:].argsort()]
            ft = np.split(faces, np.cumsum(np.unique(cols[:], return_counts=True)[1]))
            f1 = ft[0]
            f2 = ft[1]
            f3 = ft[2]
            f = [f1,f2,f3]
            poly_l = []
            zorder = [0, 1 , 0]
            alpha = [0.2,1,0.0]
            lw = [0.05,0.2,0]
            for k in range(3):
                v = [[verts_dict[ind] for ind in face] for face in f[k]]
                v = np.array(v, dtype=float)
                
                poly = Poly3DCollection(v, alpha = alpha[k],  facecolor = val[k], 
                                        edgecolors=val[k], linewidths=lw[k],zsort='max')
                
                poly_l.append(poly)    
                self.ax.add_collection3d(poly_l[-1])
                # plot plane table top
                l1 = self.lim3D[0][0]
                l2 = self.lim3D[0][1]
                X = [l1, l2, l2, l1]
                Y = [l1, l1, l2, l2]
                Z = [0,0,0,0]
                self.ax.add_collection3d(Poly3DCollection([list(zip(X, Y, Z))], zorder=0, facecolor='white'))
                #axs.set_xlim(self.lim3D[0])
                if not self.zoom:
                    # plot plane table top
                    self.ax2.add_collection3d(Poly3DCollection([list(zip(X, Y, Z))], zorder=0, facecolor='white'))
                    self.ax3.add_collection3d(Poly3DCollection([list(zip(X, Y, Z))], zorder=0, facecolor='white'))
    
                    poly2 = Poly3DCollection(v, alpha = alpha[k],  facecolor = val[k], 
                                            edgecolors=val[k], linewidths=lw[k],zsort='max')
                    self.ax2.add_collection3d(poly2)
                    poly3 = Poly3DCollection(v, alpha = alpha[k],  facecolor = val[k], 
                                        edgecolors=val[k], linewidths=lw[k],zsort='max')
                    self.ax3.add_collection3d(poly3)
                
            
            self.ax.quiver(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], color='gray')
            if not self.zoom:
                self.ax2.quiver(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], color='gray')
                self.ax3.quiver(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], color='gray')
            
            self.ax.quiver(0,0,0,v1[0], v1[1], v1[2], color='b',arrow_length_ratio=0.1)
            if not self.zoom:
                self.ax2.quiver(0,0,0,v1[0], v1[1], v1[2], color='b',arrow_length_ratio=0.1)
                self.ax3.quiver(0,0,0,v1[0], v1[1], v1[2], color='b',arrow_length_ratio=0.1)
            
            b1 = (R@v1).copy()
            self.ax.quiver(P[0], P[1], P[2], b1[0], b1[1], b1[2], color='g', linewidth=2.5, zorder=1,arrow_length_ratio=0.1)
            if not self.zoom:
                self.ax2.quiver(P[0], P[1], P[2], b1[0], b1[1], b1[2], color='g', linewidth=2.5, zorder=0,arrow_length_ratio=0.1 )
                self.ax3.quiver(P[0], P[1], P[2], b1[0], b1[1], b1[2], color='g', linewidth=2.5, zorder=0,arrow_length_ratio=0.1 )
                
            self.ax2.view_init(90, -90)
            self.ax3.view_init(0, 0)
            
        
        self.ax.quiver(0,0,0,P[0], P[1], P[2], color='r',arrow_length_ratio=0.1)
        if not self.zoom:
            self.ax2.quiver(0,0,0,P[0], P[1], P[2], color='r',arrow_length_ratio=0.1)
            self.ax3.quiver(0,0,0,P[0], P[1], P[2], color='r',arrow_length_ratio=0.1)
        
        for axs in [self.ax, self.ax2, self.ax3]:
            
            # Plot tripod
            
            tripod_mesh = mesh.Mesh.from_file(r'../data/Printed_00.stl')
            tripod_mesh.vectors = tripod_mesh.vectors
            tripod_mesh.rotate(np.array([0, 0, 1]), np.deg2rad(180))
            t = time.time()
            M_44 = np.zeros((4,4))
            M_44[:-1,-1] = P/1000
            M_44[:-1,:-1] = R
            M_44[-1,-1] = 1

            tripod_mesh.transform(M_44)
            if self.zoom:
                if axs is self.ax: # ZOOOM IN THE TOP AND SIDE VIEW
                    poly_collection = mplot3d.art3d.Poly3DCollection(tripod_mesh.vectors*1e3, zorder=10)
                    poly_collection.set_color((0.85,0.0,0.0))  # play with color
                    axs.add_collection3d(poly_collection)
            else:
                poly_collection = mplot3d.art3d.Poly3DCollection(tripod_mesh.vectors*1e3, zorder=10)
                poly_collection.set_color((0.85,0.0,0.0))  # play with color
                axs.add_collection3d(poly_collection)
            # End tripod plot
            
            # Plot endpoints and path
            axs.scatter(End_effector[0],End_effector[1],End_effector[2], color='y',zorder = 40)
            P_arr = np.array(self.P_limits)
            if self.verbose:
                print('P SHAPE', P_arr.shape)
            axs.plot(P_arr[:,0],P_arr[:,1],P_arr[:,2], 'k-', zorder = 50)
            
            axs.set_ylabel('Y (mm)')
            if axs is not self.ax3:
                if self.zoom:
                    axs.set_xlabel('X (um)')    
                    axs.set_ylabel('Y (um)')
                else:
                    axs.set_xlabel('X (mm)')    
            if axs is not self.ax2:
                if self.zoom:
                    axs.set_zlabel('Z (um)')
                    axs.set_ylabel('Y (um)')
                else:
                    axs.set_zlabel('Z (mm)')
            axs.set_xlim(self.lim3D[0])
            axs.set_ylim(self.lim3D[1])
            axs.set_zlim(self.lim3D[2])
            if self.zoom:
                if axs is not self.ax: # ZOOOM IN THE TOP AND SIDE VIEW
                    #axs.set_box_aspect((1, 1, 1), zoom=0.8)
                    #axs.margins(x=1,y=1)
                    axs.set_xlim(self.xlim)
                    
                    axs.set_ylim(self.ylim)
                    axs.set_zlim(self.zlim)
                else:
                    ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*1e-3))
                    axs.xaxis.set_major_formatter(ticks)
                    axs.yaxis.set_major_formatter(ticks)
                    axs.zaxis.set_major_formatter(ticks)
                    
            # Reescaling to mm
            else:
                ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*1e-3))
                axs.xaxis.set_major_formatter(ticks)
                axs.yaxis.set_major_formatter(ticks)
                axs.zaxis.set_major_formatter(ticks)
            
        
        self.ax2.set_zticks([])
        self.ax3.set_xticks([])
        self.ax2.set_title('Top view',y=0.9)
        self.ax3.set_title('Side view',y=0.8)

    def create_animation(self,data_2plot,name="Animation", zoom=False, mesh=False):
        
        
        # First save the mesh to plot the droplets:
        print(f'Frames: {data_2plot.shape[0]}')        
        if mesh:
            self.OUTPUT_FLAG = True
            print('Saving Meshes')
            for i in range(data_2plot.shape[0]):
                R = data_2plot[i,6:]
                X = data_2plot[i,:6]
                self.eqsystem_forward(X, *R)
                self.save_mesh(i)
            print('Meshes saved')  
        self.OUTPUT_FLAG = False
        # data_2plot is an array [i,6] with the data to plot [P,u]
        self.zoom = zoom
        R_N = R_from_vec(data_2plot[:,3:6].T)
        P_N = data_2plot[:,:3].copy()
        data_end = R_N@self.P_end + P_N
        #if self.zoom:
        mean, d = zip(*[get_mean_change(data_end, i) for i in range(3)])
        lm = np.max(d)
        self.xlim, self.ylim, self.zlim = [(m - lm, m + lm) for m in mean]
        mean = np.array([0,0,mean[-1]/2])    
        self.lim3D = [(mean[i]-self.r_tripod*2e6,mean[i]+self.r_tripod*2e6) for i in range(3)]
        self.P_limits = []
        self.fig = plt.figure(figsize=(8,5))
        self.fig.subplots_adjust(left=0.0, bottom=0, right=1, top=1, wspace=None, hspace=-0.0)
        gs = self.fig.add_gridspec(2,2,width_ratios=[2,1])
        self.ax = self.fig.add_subplot(gs[:,0], projection='3d',computed_zorder=False)
        self.ax2 = self.fig.add_subplot(gs[0,1], projection='3d',computed_zorder=False)
        self.ax2.set_proj_type('ortho')
        self.ax3 = self.fig.add_subplot(gs[1,1], projection='3d',computed_zorder=False)
        self.ax3.set_proj_type('ortho')
        ani = FuncAnimation(self.fig, self.animate,  frames=data_2plot.shape[0], repeat=False, init_func=self.init                            
                            , interval=400,fargs=(data_2plot,))  
        #ani = animation.FuncAnimation(fig, animate, frames=len(x)-1, interval=50)
        #ani.save(f'Animation.gif', fps=10, dpi=120)#, writer='imagemagick', fps=30)a  
        writervideo = FFMpegWriter(fps=5) 
        ani.save(f'../data/vids/{name}.mp4', writer=writervideo,dpi=200)
        
    
    def save_data(self):
        
        np.savetxt(f'../data/Solutions_WS_{self.V0*1e9}uL.txt', self.data, fmt='%.6e')
    
    def my_function_star(self, args):
        return eq_system_multithread(*args)    
    
    def solve_multithread(self,):
        
        return
    def plot_3D_workspace(self):
        pts = self.data[:,:3].copy()


        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        #ax.scatter(pts[:,0],pts[:,1],pts[:,2], alpha=0.1)
        
        
        verts = pts
        hull = ConvexHull(verts)
        faces = hull.simplices


        triangles = []
        for s in faces:
            #print(s)
            sq = [
                (verts[s[0], 0], verts[s[0], 1], verts[s[0], 2]),
                (verts[s[1], 0], verts[s[1], 1], verts[s[1], 2]),
                (verts[s[2], 0], verts[s[2], 1], verts[s[2], 2])
            ]
            triangles.append(sq)

        
        new_faces = WS_3D.simplify(triangles)
        for sq in new_faces:
            print(np.array(list(sq)).shape)
            triang = np.array(list(sq))
            ax.plot(triang[:,0],triang[:,1],triang[:,2], alpha=0.0)
            f = a3.art3d.Poly3DCollection([np.array(list(sq))])
            #f.set_color(colors.rgb2hex(sp.rand(3)))
            f.set_edgecolor('k')
            f.set_alpha(0.1)
            ax.add_collection3d(f)
        
        plt.show()    
        return
    
    
    def plot_test(self):
        
       # mesh = pyvista.read(r'../data/Printed_00.stl')
       # mesh.plot()
        
        
        j = 1
        i = 1
        
        file_path = f'../data/mesh/vertex_{j}_{i}.txt'
        data = pd.read_csv(file_path, sep=' ', header=None).to_numpy()
        print(data.shape)
        idx_v = data[:,0].astype(int)
        idx_n = np.arange(len(idx_v))
        #verts = data[:,1:]*1e6
        verts = data[:,:].copy()
        verts[:,1:] = verts[:,1:]*1e6
        
        
        verts[:,1] = verts[:,1]#+v1[0]
        verts[:,2] = verts[:,2]#+v1[1]
        #verts_dict = dict(zip(verts[:, 0], verts[:, 1:]))
        
        data = np.loadtxt(f'../data/mesh/facet_{j}_{i}.txt')
        
        facets = data.copy()
        faces_orig = data[:,1:4].astype(int)
        cols = data[:,-1].astype(int)
        faces = faces_orig.copy()
        faces = faces[cols[:].argsort()]
        
        
        mesh = pv.PolyData(verts, faces)
        mesh.plot(show_edges=True, line_width=5)
        

#        tripod_mesh = mesh.Mesh.from_file(r'../data/Printed_00.stl')
#        tripod_mesh.vectors = tripod_mesh.vectors
        
        
        return
    
    def pv_plot(self, data_plot, pl=None,update_cams=False, bg_plot = False):
        if data_plot.ndim == 1:
            data_plot = data_plot.reshape(1, -1)
        save_meshes(self, data_plot)
        if pl is None:
            if bg_plot:
                pl = BackgroundPlotter(window_size=(2048, 1536), shape="1|2",
                                        lighting='none')
            else:
                pl = pv.Plotter(window_size=([2048, 1536]), shape="1|2",
                            lighting='none')
            pl.enable_anti_aliasing('fxaa')
            pl.enable_lightkit()
        
        resf = data_plot[0]
        
        w = [self.w1, self.w2, self.w3]
        file_path = '../data/Printed_00.stl'
        P = resf[0:3].copy()
        u_r = resf[3:].copy()
    
        pl.subplot(0)
        plot_robot(pl,0,w,P,u_r,file_path=file_path)
        pl.camera.azimuth = 90
        if update_cams:
            pl.camera.zoom(1.5)
            cam_pos = [(-pl.camera_position[0][0], -pl.camera_position[0][1], pl.camera_position[0][2]-500),
                       (-pl.camera_position[1][0], -pl.camera_position[1][1], pl.camera_position[1][2]-500),
                       (pl.camera_position[2][0], pl.camera_position[2][1], pl.camera_position[2][2])]
            pl.camera_position = cam_pos

        #pl.camera_position = 'yz'
        pl.subplot(1)
        plot_robot(pl,0,w,P,u_r,file_path=file_path)
        #pl.add_points(End_effector)
        pl.camera_position = 'xy'
        
        pl.camera.zoom(4)
        #line = pv.Line((0, 0, 0), (0, 0, 20))
        
        pl.add_legend_scale(number_minor_ticks=2,left_axis_visibility=False,
                            bottom_axis_visibility=False,
                            right_border_offset=60,top_border_offset=60,
                            legend_visibility=False,font_size_factor=1,
                            )
        #pl.show_grid()
        pl.subplot(2)
        plot_robot(pl,0,w,P,u_r,file_path=file_path)
        pl.camera_position = 'yz'
        pl.camera.azimuth = 180
        
        pl.camera.zoom(4)
        pl.add_legend_scale(number_minor_ticks=2,left_axis_visibility=False,
                            bottom_axis_visibility=False,
                            right_border_offset=60,top_border_offset=60,
                            legend_visibility=False,font_size_factor=1,
                            )
        cam_pos = [(pl.camera_position[0][0], pl.camera_position[0][1], pl.camera_position[0][2]+300),
                   (pl.camera_position[1][0], pl.camera_position[1][1], pl.camera_position[1][2]+300),
                   (pl.camera_position[2][0], pl.camera_position[2][1], pl.camera_position[2][2])]
        pl.camera_position = cam_pos
        
        return pl
        #_ = pl.add_mesh(mesh, color='blue',show_edges=True, opacity=0.5)
        #pl.show_grid()
        
if __name__ == "__main__":
    P_end = np.array([-380,0,630]) # ANTENNA as end effector
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
    robot = Manip(P_end=P_end, V0=DROPLET_VOLUME,R_DROPLET_MIN=r_min, R_DROPLET_MAX=r_max,                    
                        ZMAX=600, ZMIN=200, R_top=150)
    # Load Data
    robot.load_data()
    r1 = 240*1e-6
    r2 = 400*1e-6
    r3 = 440*1e-6
    D1 = robot.forward_kinematics(r1,r2,r3)
    print(D1)
    data_plot = robot.data.copy()
    
    #
    frame = -1
    data_plot = D1
    print(data_plot)
    pl = robot.pv_plot(data_plot)
    pl.show()

        
