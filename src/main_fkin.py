
import numpy as np
from  manip import Manip


def set_variables():
    global R1, R2, R3
    try:
        # Prompt user for values; if empty, keep the previous value
        new_R1 = input(f"Enter new value for R1 (current: {R1}): ").strip()
        if new_R1:
            R1 = float(new_R1)
        
        new_R2 = input(f"Enter new value for R2 (current: {R2}): ").strip()
        if new_R2:
            R2 = float(new_R2)
        
        new_R3 = input(f"Enter new value for R3 (current: {R3}): ").strip()
        if new_R3:
            R3 = float(new_R3)
        
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
    robot = Manip(P_end=P_end, V0=DROPLET_VOLUME, R_DROPLET_MIN=r_min,
                  R_DROPLET_MAX=r_max, ZMAX=600, ZMIN=200, R_top=150)
    # Load Data of the foward kinematics
    robot.load_data()
    
    # ------------------
    # Example of plotting foward kinematics
    # We will set the value sof the base radious of each fluid joint
    # (R1, R2, R3)  solve the forward kinematics and plot it using pyvista
    # ----------------
    R1 = 240
    R2 = 400
    R3 = 440
    # Solve foward kinematics and save the values to D1
    D1 = robot.forward_kinematics(R1*1e-6, R2*1e-6, R3*1e-6)
    # Plot the the robot in that position
    data_plot = D1
    pl = robot.pv_plot(data_plot, bg_plot=True)

    while True:
        command = input("Enter a command: ").strip().lower()
        
        if command == "q":
            print("Exiting the program.")
            pl.close()
            break
        elif command == "c":
            pl.clear_actors()
        elif command == "s":
            set_variables()
            D1 = robot.forward_kinematics(R1*1e-6, R2*1e-6, R3*1e-6)
            data_plot = D1
            pl.clear_actors()
            robot.pv_plot(data_plot, pl)
        else:
            print("Unknown command. Try again.")



