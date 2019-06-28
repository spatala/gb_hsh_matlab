function [] = csymm_nmax()
clear all; clc;

curr_pwd = split(pwd,'/');
top_dir = '';
for ct1=1:length(curr_pwd)
    top_dir = strcat(top_dir,curr_pwd{ct1},'/');
    if (strcmp(curr_pwd{ct1},'gb_hsh_matlab'))
        break;
    end
end
util_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(util_dir));

s1 = set_vars(); Nmax = s1.Nmax; pt_grp = s1.pt_grp;

data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'cryst_symm/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ga_s, gb_s, num_gen] = get_symmgen_mats(pt_grp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
% a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1,b1);
% num_cols = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
% nsymm = size(symm_orders,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for ct4 = 1:2*num_gen
    [v0, d0] = compute_eigen(ga_s,gb_s,ct4, symm_orders);
    col0 = (abs(imag(diag(d0)))<1e-5 & abs(real(diag(d0))-1)<1e-5);
    if any(col0)
        if (ct4 == 1)
            X0 = eye(size(v0,1), size(v0,1));
        else
            X0 = orth(v1(:,col1));
        end
        [v1, d1] = combine_XY_symms(X0, v0, col0);
        col1 = (abs(imag(diag(d1)))<1e-5 & abs(real(diag(d1))-1)<1e-5);
        if ~any(col1)
            break;
        end
    end
end

if any(col1)
    save_symm_arr(symm_orders, v1, col1, data_fname0);
end
rmpath(genpath(util_dir));
end

function [v, d] = compute_eigen(ga_s,gb_s,n_gen, symm_orders)
gs1 = ga_s{n_gen}; gs2 = gb_s{n_gen};
R1 = full(generate_csymm(gs1, gs2, symm_orders));
[v,d] = eig(R1);
end

function [v, d] = combine_XY_symms(X0, v, col)
P0 = X0*X0';
Y1 = orth(v(:,col)); Q1 = Y1*Y1';
[v, d] = eig(P0*Q1*P0);
end

function save_symm_arr(symm_orders, v, col, fname)
S = orth(v(:,col));
Nmax = max(symm_orders(:));
mat_name = [fname,'Sarr_full_symm_orders_nmax_',num2str(Nmax),'.mat'];
save(mat_name,'S','symm_orders');
end