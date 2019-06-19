clear all; clc;

addpath(genpath('../Util_functions/'));
addpath(genpath('../GB_Parameters/'));

fname = get_dir_name();

Nmax = 4;

% pt_grp = 'O';
pt_grp = 'C2';

if (strcmp(pt_grp,'O'))
    num_gen = 2;
    
    a_ax1 = [0 0 1]; a_ang1 = pi/2;
    b_ax1 = [0 0 1]; b_ang1 = 0;
    
    a_ax2 = [1 0 0]; a_ang2 = pi/2;
    b_ax2 = [1 0 0]; b_ang2 = 0;
    
    a_ax3 = [0 0 1]; a_ang3 = 0;
    b_ax3 = [0 0 1]; b_ang3 = pi/2;
    
    a_ax4 = [1 0 0]; a_ang4 = 0;
    b_ax4 = [1 0 0]; b_ang4 = pi/2;
end

if (strcmp(pt_grp,'C2'))
    num_gen = 1;
    
    a_ax1 = [0 0 1]; a_ang1 = pi;
    b_ax1 = [0 0 1]; b_ang1 = 0;
    
    a_ax2 = [0 0 1]; a_ang2 = 0;
    b_ax2 = [0 0 1]; b_ang2 = pi;
end

symm_orders = zeros(Nmax^2,2);
ct3 = 1;
for ct1=0:4
    for ct2=0:4
        
        [ct1, ct2]
        
        a_val = ct1; Na = 2*a_val; 
        b_val = ct2; Nb = 2*b_val;
        c_val = min(ct1, ct2); Nc = 2*c_val;
        
        R1 = kron(rotation(a_ax1,a_ang1,Na),eye(Nb+1));
        R1 = kron(R1, eye(Nc+1));
        
        [v,d] = eig(R1);
        col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
        if any(col)
            X0 = eye(size(R1)); P0 = X0*X0';
            Y1 = orth(v(:,col)); Q1 = Y1*Y1';
            [v, d] = eig(P0*Q1*P0);
            col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
            if any(col)
                X1 = orth(v(:,col)); P1 = X1*X1';
                R2 = kron(rotation(a_ax2,a_ang2,Na),eye(Nb+1));
                R2 = kron(R2, eye(Nc+1));
                [v,d] = eig(R2);
                col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
                if any(col)
                    Y2 = orth(v(:,col)); Q2 = Y2*Y2';
                    [v, d] = eig(P1*Q2*P1);
                    col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
                    if any(col)
                        if num_gen == 1
                            save_symm_arr(a_val, b_val, v, col, pt_grp, fname);
                            symm_orders(ct3,:) = [ct1, ct2]; ct3 = ct3 + 1;
                        elseif num_gen == 2
                            X2 = orth(v(:,col)); P2 = X2*X2';
                            R3 = kron(eye(Na+1),rotation(b_ax3,b_ang3,Nb));
                            R3 = kron(R3, eye(Nc+1));
                            [v,d] = eig(R3);
                            col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
                            if any(col)
                                Y3 = orth(v(:,col)); Q3 = Y3*Y3';
                                [v, d] = eig(P2*Q3*P2);
                                col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
                                if any(col)
                                    X3 = orth(v(:,col)); P3 = X3*X3';
                                    R4 = kron(eye(Na+1),rotation(b_ax4,b_ang4,Nb));
                                    R4 = kron(R4, eye(Nc+1));
                                    [v,d] = eig(R4);
                                    col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
                                    if any(col)
                                        Y4 = orth(v(:,col)); Q4 = Y4*Y4';
                                        [v, d] = eig(P3*Q4*P3);
                                        col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
                                        if any(col)
                                            save_symm_arr(a_val, b_val, v, col, pt_grp, fname);
                                            symm_orders(ct3,:) = [ct1, ct2]; ct3 = ct3 + 1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
%
symm_orders(ct3:end,:) = [];

mat_name = [fname,'/ptgrp_',pt_grp,'/cryst_symm/symm_a_b_',pt_grp,'.mat'];
save(mat_name,'symm_orders');


rmpath(genpath('../Util_functions/'));
rmpath(genpath('../GB_Parameters/'));