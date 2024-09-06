d_scale := 1e-3 // mm to m
save:= {
        printf "" >>> "vertex.txt";
        printf "" >>> "facet.txt";
        foreach vertices do printf "%d\t %g\t %g\t %0g\n", id, x*d_scale, y*d_scale ,z*d_scale >> vertex_fname;
        foreach facet ff do 
            { printf "%d %d %d %d %d\n ",ff.id, ff.vertex[1].id,ff.vertex[2].id,ff.vertex[3].id, ff.color
                 >> facet_fname}


        }
