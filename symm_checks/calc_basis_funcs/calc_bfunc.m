%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code illustrating the computation of MBP basis functions for a set of
% random grain boundary parameters (stored in rand_gb_rots.mat)
% 
% Uses the function mbp_basis for each rotation
% Computes the full MBP basis vector.
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
pt_grp='Oh'; Nmax = 16;
coeffs_typ = 'aPLUSb_max';
% coeffs_typ = 'nmax';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Set Directory paths, add Util_functions folder

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Get data files
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,coeffs_typ,'_',num2str(Nmax),'/'];

%%%%% Load the possible (a,b) values for (pt_grp, Nmax)
mat_name = [data_fname0,'symm_ab_',pt_grp,'_',coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
a_val = symm_orders(:,1)'; b_val = symm_orders(:,2)'; c_val = min(a_val, b_val);
num_cols = sum((2*a_val+1).*(2*b_val+1).*(2*c_val+1));
%%%% tot_inds gives a mapping between (a,b) and range of indices in the
%%%% list of basis functions.
%%%% This should be made faster!
tot_inds = mbp_inds_ab_array(symm_orders);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load random GB parameters
s1 = load([top_dir,'data_files/GB_Parameters/rand_gb_rots.mat']); 
rots1 = s1.rot_mats;
n_rots = size(rots1,3);
% rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
% g1 = rots1(:,1:3); g2 = rots1(:,4:6);


%%%% Load mat file containing the symmetrized basis functions
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_gbnull_',coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
nsymm_evs = size(S,2);


ma2 = zeros(n_rots, 5);
for ct2=1:n_rots
    g2_1 = rots1(:,1:3,ct2); g2_2 = rots1(:,4:6,ct2);
    ma2(ct2,:) = rots_to_angs(g2_1, g2_2);
end

SMvec2 = zeros(n_rots,nsymm_evs);

for ct2=1:n_rots
    ct2
    Mvec2 = zeros(1,num_cols);
    for ct1 = 1:size(symm_orders,1)
        a=a_val(ct1); b = b_val(ct1);
        M2 = mbp_basis(a, b, [ma2(ct2,1), ma2(ct2,2), ...
            ma2(ct2,3), ma2(ct2,4), ma2(ct2,5)]);
        cond1 = tot_inds(:,3)==a & tot_inds(:,4)==b;
        % ind_start = find(cond1,1); ind_stop  = find(cond1,1,'last');
        tinds1 = find(cond1); ind_start = tinds1(1); ind_stop = tinds1(end);
        display([ind_start, ind_stop])
        Mvec2(ind_start:ind_stop) = M2;
    end
    SMvec2(ct2,:) = Mvec2*S;
end
