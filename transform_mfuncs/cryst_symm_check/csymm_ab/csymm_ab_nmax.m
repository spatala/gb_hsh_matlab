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

mat_name = [data_fname0,'symm_ab_',...
    pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
symm_orders = s1.symm_orders;
nsymm = size(symm_orders,1);


num_cols = 0; num_rows = 0;
for tct1=1:size(symm_orders,1)
    a_val = symm_orders(tct1,1);
    b_val = symm_orders(tct1,2);
    c_val = min(a_val, b_val);
    
    mat_name = [data_fname0,'Sarr_',...
        num2str(a_val),'_',num2str(b_val),'.mat'];
    s1 = load(mat_name);
    Svec = s1.S;
    num_cols = num_cols + size(Svec,2);
    num_rows = num_rows + (2*a_val+1)*(2*b_val+1)*(2*c_val+1);
end

S = zeros(num_rows, num_cols);

start_row_ind = 1; start_col_ind = 1;
for tct1=1:size(symm_orders,1)
    a_val = symm_orders(tct1,1);
    b_val = symm_orders(tct1,2);
    c_val = min(a_val, b_val);
    
    stop_row_ind = start_row_ind + (2*a_val+1)*(2*b_val+1)*(2*c_val+1)-1;
    
    mat_name = [data_fname0,'Sarr_',...
        num2str(a_val),'_',num2str(b_val),'.mat'];
    s1 = load(mat_name); Svec = s1.S;
    stop_col_ind = start_col_ind + size(Svec,2) - 1;
    
    S(start_row_ind:stop_row_ind, start_col_ind:stop_col_ind) = Svec;
    start_row_ind = stop_row_ind + 1;
    start_col_ind = start_col_ind + 1;
    
end

mat_name = [data_fname0,'Sarr_Nmax_',num2str(Nmax),'.mat'];
save(mat_name,'S');

rmpath(genpath(util_dir));