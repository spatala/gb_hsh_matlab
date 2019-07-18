clear all; clc;

pt_grp = 'Oh'; Nmax = 10;

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

% s1 = set_vars(); Nmax = s1.Nmax; pt_grp = s1.pt_grp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
S = s1.S;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'Sarr_MabInds_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
tot_Uprops=s1.tot_Uprops; ind_ranges=s1.ind_ranges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;
rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand()));

nrots = 1000;

rot_mats1 = rot_mats(:,:,1:nrots);
s = convert_gbrots(rot_mats1);

a1_rots  = s.a1;  a2_rots  = s.a2;
lb1_rots = s.lb1; lb2_rots = s.lb2;
q1_rots  = s.q1;  q2_rots  = s.q2;
Q1_rots  = s.Q1;  Q2_rots  = s.Q2;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Move this to a separate code.
num_rows = size(S,1);
tot_inds = mbp_inds_ab_array(symm_orders, num_rows);

% for ct1 = 1:size(S,2)
%     ct1
ct1 = 100;
st1 = ind_ranges(ct1,1);
st2 = ind_ranges(ct1,2);

U_props = tot_Uprops(st1:st2,:);
vec_inds = U_props(:,1);
a_val   = U_props(:,2); alp_val = U_props(:,5); al_val  = U_props(:,4);
b_val   =  U_props(:,3); bep_val =  U_props(:,6); be_val  = -U_props(:,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mvec = sparse(nrots,num_rows);
tic;
for ct2 = 1:size(vec_inds,1)
    U1 = rotation_wo_svd(a1_rots,lb1_rots, q1_rots, Q1_rots, a_val(ct2), alp_val(ct2), al_val(ct2));
    U2 = rotation_wo_svd(a2_rots,lb2_rots, q2_rots, Q2_rots, b_val(ct2), bep_val(ct2), be_val(ct2));
    Mvec(:,vec_inds(ct2)) = U1.*U2;
end
toc;
% end



rmpath(genpath(util_dir));