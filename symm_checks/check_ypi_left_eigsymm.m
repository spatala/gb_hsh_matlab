function [] = check_ypi_left_eigsymm()

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
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
data_fname1 = [data_fname0,'Sarr_ab/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_',coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
nsymm = size(symm_orders,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% diff1 = zeros(a_max^2,1);
% diff2 = zeros(a_max^2,1);
% diff3 = zeros(a_max^2,1);
% tct1 = 1;
s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;
rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
g1 = rots1(:,1:3); g2 = rots1(:,4:6);
ypi = vrrotvec2mat([0,1,0,pi]); g1p = ypi*g1; g2p = ypi*g2;

a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1,b1);
num_bfunc = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
Mvec_full = zeros(1, num_bfunc);
SMvec_full = zeros(1, num_bfunc);
ind_start = 1;
for ns_ord = 1:nsymm
    a_val = symm_orders(ns_ord,1); b_val = symm_orders(ns_ord,2);
    
    mat_name = [data_fname1,...
        'Sarr_ab_',num2str(a_val),'_',num2str(b_val),'_5.mat'];

    s1 = load(mat_name); S = s1.S;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Mvec = calc_Mvec(g1, g2, [a_val, b_val]);
    ind_stop = ind_start + size(Mvec,2) - 1;
    Mvec_full(1,ind_start:ind_stop) = Mvec;
    SMvec = calc_Mvec(g1p, g2p, [a_val, b_val]);
    SMvec_full(1,ind_start:ind_stop) = SMvec;
    ind_start = ind_stop + 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    norm(full(Mvec*S - SMvec*S))
    
%     diff2(tct1) = norm(Mvec*Mrot  - SMvec);
%     diff3(tct1) = norm(Mvec*Mrot1 - SMvec);
%     tct1 = tct1 + 1;
end

%     mat_name = [data_fname0,'Sarr_abc_combined_csymm_nmax_',num2str(Nmax),'.mat'];
    mat_name = [data_fname0,'Sarr_cryst_ges_nmax_',num2str(Nmax),'.mat'];
    s1 = load(mat_name); S = s1.S;
    norm(Mvec_full*S - SMvec_full*S)
end


