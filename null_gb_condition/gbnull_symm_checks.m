clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

pt_grp = 'C1'; Nmax = 1;

data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];

mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];

s1 = load(mat_name);
symm_orders = s1.symm_orders;
tot_inds = mbp_inds_ab_array(symm_orders);

mat_name = [data_fname0, ...
    'Sarr_cryst_ges_gbnull_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
S = s1.S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = load([top_dir,'data_files/GB_Parameters/rand_gb_rots.mat']); 
rot_mats = s1.rot_mats;
rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
g1 = rots1(:,1:3); g2 = rots1(:,4:6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [symm_rots, Laue] = get_symm_rots(g1,g2, pt_grp, data_fname);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_name = [data_fname0, ...
%     'Sarr_cryst_ges_gbnull_nmax_',num2str(Nmax),'.mat'];
% s1 = load(mat_name); S = (s1.S);
% nsymm_evs = size(S,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diff_vec = zeros(nsymm_evs,1);
% for ct1 = 1:nsymm_evs
% symm_Mvec = full(compute_symm_Mvec(symm_rots, S, ct1, data_fname0, Nmax));
% diff_vec(ct1) = abs(max(symm_Mvec)-min(symm_Mvec));
% % uniquetol(abs(symm_Mvec),1e-8)
% end

symm_rots = zeros(3,6,4);
symm_rots(:,:,1) = [g1, g1];
symm_rots(:,:,2) = [g1^(-1), g1^(-1)];
symm_rots(:,:,3) = [g2, g2];
symm_rots(:,:,4) = [g2^(-1), g2^(-1)];
for ct1=1:nsymm_evs 
    symm_Mvec = full(compute_symm_Mvec(symm_rots, S, ct1, data_fname0, Nmax));
end
rmpath(genpath(util_dir));
