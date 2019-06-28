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

s1 = set_vars();
Nmax = s1.Nmax;
% Nmax = 1;
pt_grp = s1.pt_grp;

data_fname0 = [top_dir,'data_files/ptgrp_',pt_grp,'/cryst_symm/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mat_name = [data_fname0,'Sarr_Nmax_',num2str(Nmax),'.mat'];
% save(mat_name,'S');
s1 = load(mat_name);
S1 = s1.S;

mat_name = [data_fname0,'Sarr_full_symm_orders.mat'];
s1 = load(mat_name); S = (s1.S);

rmpath(genpath(util_dir));