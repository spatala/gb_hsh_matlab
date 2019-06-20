function [] = basis_rot_Uab_gamma()
clear all; clc;

curr_pwd = split(pwd,'/');
top_dir = '';
for ct1=1:length(curr_pwd)
    top_dir = strcat(top_dir,curr_pwd{ct1},'/');
    if (strcmp(curr_pwd{ct1},'MATLAB_Codes'))
        break;
    end
end
util_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(top_dir));

s1 = set_vars();
Nmax = s1.Nmax; pt_grp = s1.pt_grp;

data_fname0 = [top_dir,'data_files/ptgrp_',pt_grp,'/cryst_symm/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a_ax, a_ang, b_ax, b_ang, num_gen] = get_symm_gens(pt_grp);
symm_orders = zeros(Nmax^2,2);
ct3 = 1;

for ct1=0:Nmax
    for ct2=0:Nmax
        a_val = ct1; Na = 2*a_val;
        b_val = ct2; Nb = 2*b_val;
        c_val = min(ct1, ct2); Nc = 2*c_val;
        [ct1, ct2]
        for ct4 = 1:2*num_gen
            [v0, d0] = compute_eigen(a_ax,a_ang,b_ax, b_ang, ct4, Na, Nb, Nc);
            col0 = (abs(imag(diag(d0)))<1e-5 & abs(real(diag(d0))-1)<1e-5);
            if any(col0)
                if (ct4 == 1) X0 = eye(size(v0,1), size(v0,1)); else X0 = orth(v1(:,col1)); end
                [v1, d1] = combine_XY_symms(X0, v0, col0);
                col1 = (abs(imag(diag(d1)))<1e-5 & abs(real(diag(d1))-1)<1e-5);
                if ~any(col1) break; end
            end
        end
        if any(col1)
            save_symm_arr(a_val, b_val, v1, col1, data_fname0);
            symm_orders(ct3,:) = [ct1, ct2]; ct3 = ct3 + 1;
        end        
    end
end

symm_orders(ct3:end,:) = [];
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
save(mat_name,'symm_orders');
rmpath(genpath(util_dir));
end

function [v, d] = compute_eigen(a_ax,a_ang,b_ax, b_ang, ct1, Na, Nb, Nc)
    a_ax1 = a_ax(ct1,:); b_ax1 = b_ax(ct1,:);
    a_ang1 = a_ang(ct1,:); b_ang1 = b_ang(ct1,:);
    R1 = kron(rotation(a_ax1,a_ang1,Na),rotation(b_ax1,b_ang1,Nb));
    R1 = kron(eye(Nc+1), R1);
    [v,d] = eig(R1);
end

function [v, d] = combine_XY_symms(X0, v, col) 
    P0 = X0*X0'; 
    Y1 = orth(v(:,col)); Q1 = Y1*Y1'; 
    [v, d] = eig(P0*Q1*P0);
end

function save_symm_arr(a_val, b_val, v, col, fname)
S = orth(v(:,col));
mat_name = [fname,'Sarr_',num2str(a_val),'_',num2str(b_val),'.mat'];
save(mat_name,'S');
end