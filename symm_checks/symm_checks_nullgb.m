function [] = symm_checks_nullgb(pt_grp, Nmax, coeffs_typ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Code to check the rotation of SO(4) Functions R^(a,b)
%%%%%%%
% clear all; clc;
% pt_grp = 'C1'; Nmax = 1;
% pt_grp = 'Oh'; Nmax = 4;

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

%%%% Load mat file containing the symmetrized basis functions
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_gbnull_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
nsymm_evs = size(S,2);



%%%% Set-up Identity misorientation
w_m  = 0; th_m = pi * rand(); ph_m = 2. * pi * rand();
% axang_m = [sin(th_m)*cos(ph_m), sin(th_m)*sin(ph_m), cos(th_m), w_m];
% gm = vrrotvec2mat(axang_m);

%%%% Set-up random boundary-plane orientations
w_b  = 8 * pi * (rand() - 1); th_b = pi / 2; ph_b = 2. * pi * rand();
% axang_b = [sin(th_b)*cos(ph_b), sin(th_b)*sin(ph_b), cos(th_b), w_b];
% gb = vrrotvec2mat(axang_b);



diff_vec = zeros(nsymm_evs,1);
for ct1=1:nsymm_evs
    Mvec = zeros(1,num_cols);
    for ct2 = 1:size(symm_orders,1)
        a = a_val(ct2); b = b_val(ct2);
        M1 = mbp_basis(a, b, [w_m, th_m, ph_m, w_b, ph_b]);
        
        cond1 = tot_inds(:,3)==a & tot_inds(:,4)==b;
        ind_start = find(cond1,1);
        ind_stop  = find(cond1,1,'last');
        
        Mvec(ind_start:ind_stop) = M1;

        diff_vec(ct1) = norm(Mvec*S(:,ct1));
    end
end

disp(max(abs(diff_vec)));

rmpath(genpath(util_dir));
end