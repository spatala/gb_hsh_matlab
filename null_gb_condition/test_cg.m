clear all; clc;

max_j = 2;
max_j1 = 2;
max_j2 = 2;

ct1 = 0;
ct2 = 0;
for j = 0:max_j
    for j1 = 0:max_j1
        for j2 = 0:max_j2
            for m = -j:j
                for m1 =-j1:j1
                    for m2 = -j2:j2
                        
                        c1 = clebschgordan(j1,m1,j2,m2,j,m);
                        if abs(c1) > 0
                            c1
                            [j1,m1,j2,m2,j,m]
                            ct2 = ct2 + 1;
                        end
                        ct1 = ct1 + 1;
                    end
                end
            end
        end 
    end
end
