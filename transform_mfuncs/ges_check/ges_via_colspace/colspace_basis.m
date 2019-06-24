%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code to compute the column-space basis for the GES mat.

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
addpath(genpath(top_dir));

s1 = set_vars();
Nmax = s1.Nmax; pt_grp = s1.pt_grp;

data_fname0 = [top_dir,'data_files/ptgrp_',pt_grp,'/cryst_symm/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

symm_orders = zeros(Nmax^2,2);
ct3 = 1;
for ct1=0:Nmax
    for ct2=0:Nmax
        symm_orders(ct3,:) = [ct1,ct2];
        ct3 = ct3 + 1;
    end
end

a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1,b1);
num_rows = sum((2*a1+1).*(2*b1+1).*(2*c1+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% The variable tot_inds contains the mapping between the column-index
%%%% for the M-function and the (a,b,gamma,alpha, beta) values.
tot_inds = mbp_inds_ab_array(symm_orders, num_rows);

%%%% ges_mat is the $A$ matrix that specifies the grain-exchange-symmetry
%%%% for the coefficients.
num_cols = num_rows;
ges_mat = sparse(num_rows, num_cols);
for ct1=1:num_cols
%%%% Get the (a,b,gamma,alpha, beta) corresponding to the column.
    a1 = tot_inds(ct1,3); b1 = tot_inds(ct1,4);
    gamma1 = tot_inds(ct1,5); 
    alpha1 = tot_inds(ct1,6); beta1 = tot_inds(ct1,7);
    
%%%% Find the index for (b,a,gamma,beta, alpha) row.
    a2 = b1; b2 = a1;
    gamma2 = gamma1;
    alpha2 = beta1; beta2 = alpha1;
    ind1 = find((tot_inds(:,3) == a2) & ...
        (tot_inds(:,4) == b2) & ...
        (tot_inds(:,5) == gamma2) & ...
        (tot_inds(:,6) == alpha2) & ...
        (tot_inds(:,7) == beta2));
    
%%%% The (a,b, gamma, alpha, beta)th column contains a 1 in the 
%%%% (b,a, gamma, beta, alpha)th row
    ges_mat(ind1,ct1) = 1;
%%%% The (a,b, gamma, alpha, beta)th column contains a (?1)^{a+b} in the 
%%%% (a,b, gamma, alpha, beta)th row
    ges_mat(ct1,ct1) = (-1)^(a1+b1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Finding the column-space for a sparse matrix
col1 = get_colspace(ges_mat);
%%%% Save the column-space in .mat file.
mat_name = ['Y_ges_Nmax_',num2str(Nmax),'.mat']; save(mat_name, 'col1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coeffs = col1*(rand(1749,1)*100);

s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;
rots1 = rot_mats(:,:,1); r1 = rots1(:,1:3); r2 = rots1(:,4:6);
M1 = calc_Mfunc_ab_array(symm_orders, num_rows, rots1);
nrots = size(rot_mats,3);

g_ypi = vrrotvec2mat([0,1,0, pi]);
r1_ges = g_ypi*r2; r2_ges = g_ypi*r1;
rots_ges = [r1_ges, r2_ges];
M1_ges = calc_Mfunc_ab_array(symm_orders, num_rows, rots_ges);

