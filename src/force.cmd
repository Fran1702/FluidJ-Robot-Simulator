// force.cmd

// Force calculations.  Using the principle of virtual work,
//  Central differences are
// used for greater accuracy.  It is assumed that forces are calculated
// at an equilibrium, so the surface does not have to be re-evolved
// after each perturbation.  Instead, the pressure can be used to 
// compensate for the volume change, since pressure is the Lagrange
// multiplier for the volume constraint.  The surface is left in its
// original position after each calculation.
// Note: To check these force calculations, you can do a perturbation
// by hand with evolution.

// This file should be read in after change_simple.cmd

f_scale:= 1e-12 // Note Energy in Joule-> 1e12 in mg*mm^2/s^2 and displacement already in mm
t_scale:= 1e-12

calc_xf := { f_dx := RMAX/40000;  // small shift
             new_x0 := x0 + f_dx; change_x0; 
             energy_hi := total_energy - body[1].pressure*(body[1].volume-body[1].target);
             new_x0 := x0 - 2*f_dx; change_x0; 
             energy_lo := total_energy - body[1].pressure*(body[1].volume-body[1].target);
             xforce := -(energy_hi - energy_lo)/(2*f_dx)*f_scale;
             new_x0 := x0 + f_dx; change_x0; 
             printf "xforce: %17.15g\n",xforce;
           }

calc_yf := { f_dy := RMAX/40000;  // small shift
             new_y0 := y0 + f_dy; change_y0; 
             energy_hi := total_energy - body[1].pressure*(body[1].volume-body[1].target);
             new_y0 := y0 - 2*f_dy; change_y0; 
             energy_lo := total_energy - body[1].pressure*(body[1].volume-body[1].target);
             yforce := -(energy_hi - energy_lo)/(2*f_dy)*f_scale;
             new_y0 := y0 + f_dy; change_y0; 
             printf "yforce: %17.15g\n",yforce;
           }

calc_zf := { f_dz := RMAX/40000;  // small shift
             new_z0 := z0 + f_dz; change_z0; 
             energy_hi := total_energy - body[1].pressure*(body[1].volume-body[1].target);
             new_z0 := z0 - 2*f_dz; change_z0; 
             energy_lo := total_energy - body[1].pressure*(body[1].volume-body[1].target);
             zforce := -(energy_hi - energy_lo)/(2*f_dz)*f_scale;
             new_z0 := z0 + f_dz; change_z0; 
             printf "zforce: %17.15g\n",zforce;
           }

// Torque calculations
calc_yt := { t_dang := 1/1000;  // small shift
             new_ang_y := ang_y + t_dang; change_ang_y; 
             energy_hi := total_energy - body[1].pressure*(body[1].volume-body[1].target);
             new_ang_y := ang_y - 2*t_dang; change_ang_y; 
             energy_lo := total_energy - body[1].pressure*(body[1].volume-body[1].target);
             ytorque := -(energy_hi - energy_lo)/(2*t_dang)*t_scale;
             new_ang_y := ang_y + t_dang; change_ang_y; 
             printf "ytorque: %17.15g\n",ytorque;
           }

calc_xt := { t_dang := 1/1000;  // small shift
             new_ang_x := ang_x + t_dang; change_ang_x; 
             energy_hi := total_energy - body[1].pressure*(body[1].volume-body[1].target);
             new_ang_x := ang_x - 2*t_dang; change_ang_x; 
             energy_lo := total_energy - body[1].pressure*(body[1].volume-body[1].target);
             xtorque := -(energy_hi - energy_lo)/(2*t_dang)*t_scale;
             new_ang_x := ang_x + t_dang; change_ang_x; 
             printf "xtorque: %17.15g\n",xtorque;
           }
           
calc_zt := { t_dang := 1/1000;  // small shift
             new_ang_z := ang_z + t_dang; change_ang_z; 
             energy_hi := total_energy - body[1].pressure*(body[1].volume-body[1].target);
             new_ang_z := ang_z - 2*t_dang; change_ang_z; 
             energy_lo := total_energy - body[1].pressure*(body[1].volume-body[1].target);
             ztorque := -(energy_hi - energy_lo)/(2*t_dang)*t_scale;
             new_ang_z := ang_z + t_dang; change_ang_z; 
             printf "ztorque: %17.15g\n",ztorque;
           }

//  Handy command to calculate all the forces and torques
calc_all := { calc_xf; calc_yf; calc_zf; calc_yt; calc_xt; }

