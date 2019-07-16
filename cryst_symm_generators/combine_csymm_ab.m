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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname0 = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'nmax_',num2str(Nmax),...
    '/symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsymm = size(symm_orders,1);
n_eigv = 0;
for ct1 = 1:nsymm
    a_val = symm_orders(ct1,1); Na = 2*a_val;
    b_val = symm_orders(ct1,2); Nb = 2*b_val;
    c_val = min(a_val, b_val);  Nc = 2*c_val;
    
    mat_name = ['symm_ab_',pt_grp,'_a_', ...
                num2str(a_val), '_b_', num2str(b_val), '.mat'];
    s1 = load(mat_name);
    S = s1.S;
    S1 = kron(eye(Nc+1),S);
    n_eigv = n_eigv+size(S1,2);
end

