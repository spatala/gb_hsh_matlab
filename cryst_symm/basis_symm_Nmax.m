clear all; clc;

addpath(genpath('../Util_functions/'));

Nmax = 4;

fname = get_dir_name();

% pt_grp = 'O';
pt_grp = 'C2';
mat_name = [fname,'/ptgrp_',pt_grp,'/cryst_symm/symm_a_b_',pt_grp,'.mat'];
s1 = load(mat_name);
symm_orders = s1.symm_orders;

num_cols = 0; num_rows = 0;
for tct1=1:size(symm_orders,1)
    a_val = symm_orders(tct1,1);
    b_val = symm_orders(tct1,2);
    c_val = min(a_val, b_val);
    
    mat_name = [fname,'/ptgrp_',pt_grp,'/cryst_symm/Sarr_',num2str(a_val),'_',num2str(b_val),'.mat'];
    s1 = load(mat_name);
    Svec = s1.S;
    num_cols = num_cols + size(Svec,2);
    num_rows = num_rows + (2*a_val+1)*(2*b_val+1)*(2*c_val+1);
end

Sarr_nmax = zeros(num_rows, num_cols);

start_row_ind = 1;
start_col_ind = 1;
for tct1=1:size(symm_orders,1)
    a_val = symm_orders(tct1,1);
    b_val = symm_orders(tct1,2);
    c_val = min(a_val, b_val);
    
    stop_row_ind = start_row_ind + (2*a_val+1)*(2*b_val+1)*(2*c_val+1)-1;
    
    mat_name = [fname,'/ptgrp_',pt_grp,'/cryst_symm/Sarr_',num2str(a_val),'_',num2str(b_val),'.mat'];
    s1 = load(mat_name);
    Svec = s1.S;
    stop_col_ind = start_col_ind + size(Svec,2) - 1;
    
    Sarr_nmax(start_row_ind:stop_row_ind, start_col_ind:stop_col_ind) = Svec;
    start_row_ind = stop_row_ind + 1;
    start_col_ind = start_col_ind + 1;
    
end

mat_name = [fname,'/ptgrp_',pt_grp,'/cryst_symm/Sarr_Nmax_',num2str(Nmax),'_',pt_grp,'.mat'];
save(mat_name,'Sarr_nmax');
rmpath(genpath('../Util_functions/'));