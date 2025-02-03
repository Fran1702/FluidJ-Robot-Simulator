import numpy as np
from  manip import Manip
import matplotlib.pyplot as plt

def set_variables():
    global R1, R2, R3
    try:
        # Prompt user for values; if empty, keep the previous value
        new_R1 = input(f"Enter new value for R1 (current: {R1}, range: {r_min}-{r_max}): ").strip()
        if new_R1:
            R1_candidate = float(new_R1)
            if r_min <= R1_candidate <= r_max:
                R1 = R1_candidate
            else:
                print(f"Value out of range! R1 must be between {r_min} and {r_max}.")
        
        new_R2 = input(f"Enter new value for R2 (current: {R2}, range: {r_min}-{r_max}): ").strip()
        if new_R2:
            R2_candidate = float(new_R2)
            if r_min <= R2_candidate <= r_max:
                R2 = R2_candidate
            else:
                print(f"Value out of range! R2 must be between {r_min} and {r_max}.")
        
        new_R3 = input(f"Enter new value for R3 (current: {R3}, range: {r_min}-{r_max}): ").strip()
        if new_R3:
            R3_candidate = float(new_R3)
            if r_min <= R3_candidate <= r_max:
                R3 = R3_candidate
            else:
                print(f"Value out of range! R3 must be between {r_min} and {r_max}.")
        
        print(f"Updated variables: R1={R1}, R2={R2}, R3={R3}")
    except ValueError:
        print("Invalid input! Please enter valid numerical values.")



if __name__ == "__main__":
    # ------------------
    # Robot definitions
    # ------------------
    P_end = np.array([-380, 0, 630]) # ANTENNA as end effector
    DROPLET_VOLUME = 0.102e-9        # Volume of the fluid joints
    theta_0 = 135*np.pi/180          # Max contact angle (to compute r_min)
    theta_min = 60*np.pi/180         # Min contact angle (to compute r_max)
    # Computation of the max and min values of the base radius 
    # of the flui joints
    r_min = (3*DROPLET_VOLUME*np.sin(theta_0)**2 /
             (np.pi*(2+np.cos(theta_0))*(1-np.cos(theta_0))**2))**(1/3)*1e6
    r_max = (3*DROPLET_VOLUME*np.sin(theta_min)**2 /
             (np.pi*(2+np.cos(theta_min))*(1-np.cos(theta_min))**2))**(1/3)*1e6
    # Round
    r_min = np.round(r_min*1.0, 0)
    r_max = np.round(r_max*1.052, 0)
    # create the mani object (robot) with the desired settings
    # Zmax and Zmin are sets as max and min values for the z displacement
    # And used only to avoid going to far
    w1 = np.array([800,0,0])
    w2 = np.array([0,400,0])
    w3 = np.array([-800,0,0])
    w = [w1,w2,w3]
    robot = Manip(p_end=P_end, vol=DROPLET_VOLUME, r_droplet_min=r_min,
                  w=w, z_max=600, z_min=200, r_top=150,
                  stl_name = 'leg')
                  #r_droplet_max=r_max, z_max=600, z_min=200, r_top=150)
    # Load Data of the foward kinematics
    robot.load_data()
    
    # ------------------
    # Example of plotting foward kinematics
    # We will set the value sof the base radious of each fluid joint
    # (R1, R2, R3)  solve the forward kinematics and plot it using pyvista
    # ----------------
    R1 = 350
    R2 = 350
    R3 = 350
    # Solve foward kinematics and save the values to D1
    D1 = robot.forward_kinematics(R1*1e-6, R2*1e-6, R3*1e-6)
    # Plot the the robot in that position
    data_plot = D1
    #R_arr = 100*np.sin(np.linspace(0,2*np.pi,5))+R2
    delta = 50*np.sin(np.linspace(0.5*np.pi,1.5*np.pi,30))
    K = 0.4
    r1 = np.ones_like(delta)*R2 - delta + K*np.abs(delta)
    r2 = np.ones_like(delta)*R2 + delta + K*np.abs(delta)
    r3 = np.ones_like(delta)*R2 - delta + K*np.abs(delta)
    # swing phase
    d2 = 50*np.sin(np.linspace(1.5*np.pi,2.5*np.pi,10)) 
    sr1 = np.ones_like(d2)*R2 - d2 - K*np.abs(d2) - 20
    sr2 = np.ones_like(d2)*R2 + d2 - K*np.abs(d2) - 20
    sr3 = np.ones_like(d2)*R2 - d2 - K*np.abs(d2) - 20
    rt1 = np.concatenate((r1,sr1))
    rt2 = np.concatenate((r2,sr2))
    rt3 = np.concatenate((r3,sr3))
#    plt.plot(rt1)
#    plt.plot(rt2)
#    plt.plot(rt3)
#    plt.show()
#    input()
    pl = robot.pv_plot(data_plot, bg_plot=True, save_stl=True)
    for i in range(len(rt1)):
            print(f'i: {i}, Radius: {rt1[i]}, {rt2[i]}, {rt3[i]}')
            D1 = robot.forward_kinematics(rt1[i]*1e-6,
                                          rt2[i]*1e-6,
                                          rt3[i]*1e-6)
            data_plot = D1
            pl.clear_actors()
            robot.pv_plot(data_plot, pl, save_stl=True)





