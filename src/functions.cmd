
// Typical evolution
gogo := { u; r; g 5; r; g 5; r; g 5; hessian; hessian }


// A farther evolution, to 
gofar := { gogo; r; g 10 ; hessian; hessian; r; g 10; hessian; hessian }
      
//check_increase ON

re := {refine edges where on_constraint 1 and length>lmin }

gogo := {
            { 
            refine edge where length > 2*lmin and not no_refine;
            g 5;
            V 3; u; V 5;
            refine edge where length > 2*lmin and not no_refine;
            g 10; 
            V 3; u; V 3;
            U; g 10;
            hessian_seek; hessian_seek;
            };
          }

script := { 
            gogo;
            //t 0.95*lmin; 
            refine edge where length > 2*lmin and not no_refine;
            //r;
            //u; V 5; u; V 5; u; g 3;
            u; V 2; u; V 2; g 3;
            hessian_seek;
            {
            
            refine edge where length > 2*lmin and not no_refine;
	    u;
            t 0.95*lmin; 
            //u; V 4; u; V 3; g 5; hessian; hessian; 
            u; V 10; g 10; hessian_seek;hessian_seek;  
            } 2;
            u;V 1; g 3;// hessian;hessian;
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
