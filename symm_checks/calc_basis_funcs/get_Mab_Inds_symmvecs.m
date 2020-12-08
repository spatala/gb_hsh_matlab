clear all; clc;

pt_grp = 'Oh'; Nmax = 16; coeffs_typ = 'aPLUSb_max';

top_dir = get_top_dir();
util_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(util_dir));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,coeffs_typ,'_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_',coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_gbnull_',coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
S = s1.S;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'Sarr_MabInds_',coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
tot_Uprops=s1.tot_Uprops; ind_ranges=s1.ind_ranges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = load([top_dir,'data_files/GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;
rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand()));

nrots = 10;

rot_mats1 = rot_mats(:,:,1:nrots);
s = convert_gbrots(rot_mats1);

a1_rots  = s.a1;  a2_rots  = s.a2;
lb1_rots = s.lb1; lb2_rots = s.lb2;
q1_rots  = s.q1;  q2_rots  = s.q2;
Q1_rots  = s.Q1;  Q2_rots  = s.Q2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Move this to a separate code.
num_rows = size(S,1);
tot_inds = mbp_inds_ab_array(symm_orders);

% for ct1 = 1:size(S,2)
%     ct1
%     % ct1 = 8;
%     st1 = ind_ranges(ct1,1);
%     st2 = ind_ranges(ct1,2);
%     
%     U_props = tot_Uprops(st1:st2,:);
%     vec_inds = U_props(:,1);
%     a_val   = U_props(:,2); alp_val = U_props(:,5); al_val  = U_props(:,4);
%     b_val   =  U_props(:,3); bep_val =  U_props(:,6); be_val  = -U_props(:,4);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     Mvec = sparse(nrots,num_rows);
%     tic;
%     for ct2 = 1:size(vec_inds,1)
%         U1 = rotation_wo_svd(a1_rots,lb1_rots, q1_rots, Q1_rots, a_val(ct2), alp_val(ct2), al_val(ct2));
%         U2 = rotation_wo_svd(a2_rots,lb2_rots, q2_rots, Q2_rots, b_val(ct2), bep_val(ct2), be_val(ct2));
%         Mvec(:,vec_inds(ct2)) = U1.*U2;
%     end
%     toc;
% end


Mvec1 = compute_Mvec(s, nrots, num_cols, vec_inds);


rmpath(genpath(util_dir));