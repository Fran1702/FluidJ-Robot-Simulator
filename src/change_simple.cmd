// change_simple.cmd

new_x0 := x0
change_x0 := { change_dx := new_x0 - x0; x0 := new_x0;
               set vertex x x+change_dx*1e3 where fixed and z > z0*1e3/2; // TOP SURFACE
               recalc;
             }

new_y0 := y0
change_y0 := { change_dy := new_y0 - y0; y0 := new_y0;
               set vertex y y+change_dy*1e3 where fixed and z > z0*1e3/2;
               recalc;
             }
new_z0 := z0
change_z0 := { change_dz := new_z0 - z0; old_z0 := z0; z0 := new_z0;
               set vertex z z+change_dz*1e3 where fixed and z > z0*1e3/2 ;
               recalc;
             }

new_ang_y := ang_y
change_ang_y := { change_angy := new_ang_y - ang_y; old_ang_y := ang_y; ang_y := new_ang_y;
               recalc;
             }


new_ang_x := ang_x
change_ang_x := { change_angx := new_ang_x - ang_x; old_ang_x := ang_x; ang_x := new_ang_x;
               recalc;
             }

new_ang_z := ang_z
change_ang_z := { change_angz := new_ang_z - ang_z; old_ang_z := ang_z; ang_z := new_ang_z;
               recalc;
             }
