
// Typical evolution
gogo := { u; r; g 5; r; g 5; r; g 5; hessian; hessian }


// A farther evolution, to 
gofar := { gogo; r; g 10 ; hessian; hessian; r; g 10; hessian; hessian }
      
//check_increase ON

re := {refine edges where on_constraint 1 and length>lmin }


//gg := { g; 
//	change_energy := (total_energy-olde)/total_energy;
	//printf "Change: %10.5g\n",change_energy; olde:=total_energy}

//gg2:= {while abs(change_energy)>0.1 do {g;V;u;V;olde := total_energy;gg;};}
olde:=total_energy;

pr :={printf "ENE: %10.15e \n", total_energy;}

actual_mean_l := 1;
get_mean_l := {
	val := 0;
	k_n := 0;
	foreach edge where not no_refine do {val:= val + length; k_n :=k_n+1; };
	mean := val/k_n;
//	printf "mean: %1.5g \n", mean;
	actual_mean_l := mean;
	}
hlength :={histogram(edges where not no_refine,length):}

gogo := {
	t0 := clock;
            {
            refine edge where length > 2*lmin and not no_refine;
            g; V 5;u; V 5; g 5;
            it_ref := 0;
            quadratic;
	     while edge_refine_count/edge_count || edge_delete_count/edge_count > 0.1 do
            {	
//	    printf "Deleted: %10.5g \n", edge_delete_count;
//	    printf "refine: %10.5g \n", edge_refine_count;
	         printf "Refination: %10.5g\n", it_ref; 
	         it_ref := it_ref+1;
 	         if it_ref == 2 then conj_grad on;
	         it_g :=0;
	            change_energy_TENS :=1;
	            while abs(change_energy_TENS) > 1e-4 do {
            	        olde := total_area;
	            	it_g :=it_g+1;
	            	printf "Iteration: %10.5g\n", it_g;
	            	printf "Elapsed time before g: %10.5g seconds\n", clock-t0;
	           	u;u; g 1 ;
	            	if it_g%5 == 0 then {hessian_seek;};
	           	//printf "Elapsed time after g: %10.5g seconds\n", clock-t0;
            		// printf "Change: %10.5g\n",(total_energy-olde)/total_energy;
            		change_energy_TENS := (total_area-olde)/total_area;
//                	//printf "ChangeEnergy/TENS: %10.5e \n", abs(change_energy_TENS);

	            	if it_g == 50 then break;
	            	
	            	};
 	           if it_ref == 8 then break;
	
        	  reset_counts;
      	      //   printf "Elapsed time before refine: %10.5g seconds\n", clock-t0; 
	      get_mean_l; // computes the mean of the legth (edges)
     	      if actual_mean_l/2 > lmin then {t 0.95*actual_mean_l/2;};
	      if actual_mean_l/2 < 0.9*lmin then  t 0.9*lmin;
              refine edge where length > 2*lmin and not no_refine;
//     	       printf "Elapsed time after refine: %10.5g seconds\n", clock-t0;
      	      //   printf "Elapsed time before pop edges: %10.5g seconds\n", clock-t0;
     	      //   printf "Elapsed time after pop edges: %10.5g seconds\n", clock-t0;
            };
          };
          printf "Elapsed before hessians: %10.5g seconds\n", clock-t0;
//         e;
	 olde := total_energy;
	 printf "Energy: %10.5g \n", total_energy;
	  hessian;//_seek; hessian;
	  printf "Energy: %10.5g \n", total_energy;
	  printf "Change: %10.5g\n",(total_energy-olde)/total_energy;
//	  e;
          printf "Elapsed time: %10.5g seconds\n", clock-t0;
          }

script := { 
            gogo;
           // refine edge where length > 2*lmin and not no_refine;
           // u; V 2; u; V 2; g 3;
           // hessian_seek;
            // {
           // refine edge where length > 2*lmin and not no_refine;
	    // u;
           // t 0.95*lmin; 
           // u; V 10; g 10; hessian_seek;hessian_seek;  
           // } 2;
            //u;V 1; g 3;// hessian;hessian;
          }
t_dene :={gogo;
	printf "Initial energy: ";
	pr;
	init_ene:= total_energy;
	new_x0 := 1e-9; // Small displacement
	change_x0; 
	 gogo;
 	printf "Energy after displacement: ";
 	pr;
	printf "Energy change: %1.10e \n ", total_energy-init_ene;
		}

Calc_area := {base_area := sum(facet where color==white,area);};

run_forces := {
                script;
                printf "\n";
                calc_xf;
                calc_yf;
                calc_zf;
                calc_xt;
                calc_yt;
                calc_zt;
                }
                

