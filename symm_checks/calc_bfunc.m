%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculate MBP basis function (a,b <= Nmax)
%%%% 
%%%% 
%%%% 
%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
pt_grp='Oh'; Nmax = 16;
coeffs_typ = 'aPLUSb_max';

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
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];

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
%%%% Load random GB parameter
s1 = load([top_dir,'data_files/GB_Parameters/rand_gb_rots.mat']); 
symm_rots = s1.rot_mats;
nsymm_rots = size(symm_rots,3);
% rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
% g1 = rots1(:,1:3); g2 = rots1(:,4:6);


%%%% Load mat file containing the symmetrized basis functions
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_gbnull_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
nsymm_evs = size(S,2);

%%%% Compute the basis functions for all the symmetrically equivalent
%%%% representations of the random GB parameter.
ma2 = zeros(nsymm_rots, 5);
for ct2=1:nsymm_rots
    g2_1 = symm_rots(:,1:3,ct2); g2_2 = symm_rots(:,4:6,ct2);
    ma2(ct2,:) = rots_to_angs(g2_1, g2_2);
end

SMvec2 = zeros(nsymm_rots,nsymm_evs);

for ct2=1:nsymm_rots
    ct2
    Mvec2 = zeros(1,num_cols);
    for ct1 = 1:size(symm_orders,1)
        a=a_val(ct1); b = b_val(ct1);
        M2 = mbp_basis(a, b, [ma2(ct2,1), ma2(ct2,2), ...
            ma2(ct2,3), ma2(ct2,4), ma2(ct2,5)]);
        cond1 = tot_inds(:,3)==a & tot_inds(:,4)==b;
        ind_start = find(cond1,1); ind_stop  = find(cond1,1,'last');
        Mvec2(ind_start:ind_stop) = M2;
    end
    SMvec2(ct2,:) = Mvec2*S;
end
