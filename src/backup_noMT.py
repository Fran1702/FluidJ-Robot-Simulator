#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 10:54:24 2024

@author: hector.ortiz
"""

    def Jacob_inverse(self,Q, *args):
        # Computes the jacobian for the inverse kinematics
        d = 1e-5
        JM = np.zeros((6,6))
        rx, ry, rz = args
        # create the thread pool
        Q_l = np.array([Q.copy() for i in range(6)]) + np.eye(6)*d
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
        
        for i  in range(6):
            JM[:,i] = dF_l[i]/d
        return JM
    
    
    def eqsystem_direct(self, Q, *args):
        rbase1,rbase2,rbase3 = args
        if self.verbose:
            print(Q)
            
          #  input()
        #input()
        Q = Q.copy()
        # Q is {u{3},P{3}}
        # rb is the base radius of each droplet
        #u_vec = np.ones(3)
        #u_vec[0:2] = Q[0:2]
        u_vec = Q[3:]
        u_vec = u_vec - self.U_OFFSET_SOL  # scale the values to improve solution
        R = R_from_vec(u_vec)
        a_x, a_y, a_z = rot_2_angles(R, deg=False)
        P = Q[0:3]
        P[0] = P[0]-self.OFF_SET_XY_SOL
        P[1] = P[1]-self.OFF_SET_XY_SOL
        P[2] = P[2]-self.Z_OFFSET_SOL # Remove the offset
        

        
        # Kinematics equation [um]
        H1 = R@self.b1 + P - self.w1
        H2 = R@self.b2 + P - self.w2
        H3 = R@self.b3 + P - self.w3
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
        n_tripod = np.cross((R@self.b1-R@self.b2),(R@self.b3-R@self.b2)) # Normal vector to the plane of the tripod (x_t, y_x)
        
        # The following cycle can be improved using threads...
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
                if self.verbose:
                    print('D is zero')
            # Rotatation matrix  to rotate in z axis 

            R_se = R@(R_ang.T)
            
            phi, theta,psi = rot_2_angles(R_se)
            ang_x = phi 
            ang_y = theta
            ang_z = psi
            if self.verbose:
                txt = rf'φx: {ang_x:.2f}$^\circ$, θy: {ang_y:.2f}$^\circ$, ψz: {ang_z:.2f}$^\circ$'
                print(txt)
            ang_x = phi*np.pi/180
            ang_y = theta*np.pi/180
            ang_z = psi*np.pi/180

            fx, fy, fz, tx, ty, tz = calc_forces(x,y,z, a_x, a_y ,a_z, rbase, outputfiles=self.OUTPUT_FLAG)
            # projection to original xy axis
            f.append(np.array([fx,fy,fz]))
            
            t.append(np.array([tx,ty,tz]))
            #if self.verbose:
             #   input('i')
            
        # Force equation:
        F = f[0]+f[1]+f[2] # [N]
        F = F*1e6 # To uN
        if self.verbose:
            print('F: ',F)
        
        # Torque equations f in N b in um:
        T = np.cross(R@self.b1, f[0]) + np.cross(R@self.b2, f[1]) + np.cross(R@self.b3, f[2])
        T = T+t[0]*1e6+t[1]*1e6+t[2]*1e6
        if self.verbose:
            print('Torq vecs 1: ',np.cross(R@self.b1+P, f[0]))
            print('Torq vecs 2: ',np.cross(R@self.b2+P, f[1]) )
            print('Torq vecs 3: ', np.cross(R@self.b3+P, f[2]))
            print('t1,t2,t3', t[0], t[1], t[2])
            print('T: ',T)
            
        T = T # [uNm]
        return np.concatenate((F,T[:3])).flatten()
    
    
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
        R = R_from_vec(u_vec)
        a_x, a_y, a_z = rot_2_angles(R, deg=False)
        # transform P_end to the translation of the frame reference P
        P = r_end - R@self.P_end

    
        # Kinematics equation [um]
        H1 = R@self.b1 + P - self.w1
        H2 = R@self.b2 + P - self.w2
        H3 = R@self.b3 + P - self.w3
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
        n_tripod = np.cross((R@self.b1-R@self.b2),(R@self.b3-R@self.b2)) # Normal vector to the plane of the tripod (x_t, y_x)
        
        # The following cycle can be improved using threads...
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
                if self.verbose:
                    print('D is zero')
            # Rotatation matrix  to rotate in z axis 

            R_se = R@(R_ang.T)
            
            phi, theta,psi = rot_2_angles(R_se)
            ang_x = phi 
            ang_y = theta
            ang_z = psi
            if self.verbose:
                txt = rf'φx: {ang_x:.2f}$^\circ$, θy: {ang_y:.2f}$^\circ$, ψz: {ang_z:.2f}$^\circ$'
                print(txt)
            ang_x = phi*np.pi/180
            ang_y = theta*np.pi/180
            ang_z = psi*np.pi/180

            fx, fy, fz, tx, ty, tz = calc_forces(x,y,z, a_x, a_y ,a_z, rbase, outputfiles=self.OUTPUT_FLAG)
            # projection to original xy axis
            f.append(np.array([fx,fy,fz]))
            
            t.append(np.array([tx,ty,tz]))
            #if self.verbose:
             #   input('i')
            
        # Force equation:
        F = f[0]+f[1]+f[2] # [N]
        F = F*1e6 # To uN
        if self.verbose:
            print('F: ',F)
        
        # Torque equations f in N b in um:
        T = np.cross(R@self.b1, f[0]) + np.cross(R@self.b2, f[1]) + np.cross(R@self.b3, f[2])
        T = T+t[0]*1e6+t[1]*1e6+t[2]*1e6
        if self.verbose:
            print('Torq vecs 1: ',np.cross(R@self.b1+P, f[0]))
            print('Torq vecs 2: ',np.cross(R@self.b2+P, f[1]) )
            print('Torq vecs 3: ', np.cross(R@self.b3+P, f[2]))
            print('t1,t2,t3', t[0], t[1], t[2])
            print('T: ',T)
            
        T = T # [uNm]
        return np.concatenate((F,T[:3])).flatten()    