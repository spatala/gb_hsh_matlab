%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Code to check the rotation of SO(4) Functions R^(a,b)
%%%%%%%
clear all; clc;

pt_grp = 'Oh'; Nmax = 14;

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

% s1 = set_vars();
% Nmax = s1.Nmax; pt_grp = s1.pt_grp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = load([top_dir,'data_files/GB_Parameters/rand_gb_rots.mat']); 
rot_mats = s1.rot_mats;
rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
g1 = rots1(:,1:3); g2 = rots1(:,4:6);
% M1 = vrrotvec2mat([1,1,1,pi/3]);
% th = rand()*pi; ph = rand()*2*pi;
% nvec = [sin(th)*cos(ph); sin(th)*sin(ph); cos(th)];
% rots1 = mbp_to_rots([M1, nvec]);
% g1 = rots1(:,1:3); g2 = rots1(:,4:6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[symm_rots, Laue] = get_symm_rots(g1,g2, pt_grp, data_fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
nsymm_evs = size(S,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ct1 = 1:nsymm_evs
symm_Mvec = full(compute_symm_Mvec(symm_rots, S, ct1, data_fname0, Nmax));
uniquetol(abs(symm_Mvec),1e-8)
end

rmpath(genpath(util_dir));