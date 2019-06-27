%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Code to check the rotation of SO(4) Functions R^(a,b)
%%%%%%%
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

s1 = set_vars();
Nmax = s1.Nmax; pt_grp = s1.pt_grp;

data_fname0 = [top_dir,'data_files/ptgrp_',pt_grp,'/cryst_symm/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symm_orders = zeros(Nmax^2,2);
ct3 = 1;
for ct1=0:Nmax
    for ct2=0:Nmax
        symm_orders(ct3,:) = [ct1,ct2];
        ct3 = ct3 + 1;
    end
end
a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1,b1);
num_cols = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
nsymm = size(symm_orders,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;
rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
g1 = rots1(:,1:3); g2 = rots1(:,4:6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zpi = vrrotvec2mat([0,0,1,pi]);

% %%%%
% gs1 = eye(3); gs2 = eye(3);
% Sg1 = g1*gs1; Sg2 = g2*gs2;

% %%%%
% gs1 = Zpi; gs2 = eye(3);
% Sg1 = g1*gs1; Sg2 = g2*gs2;

% %%%%
% gs1 = eye(3); gs2 = Zpi;
% Sg1 = g1*gs1; Sg2 = g2*gs2;
% 

% %%%%
% gs1 = Zpi; gs2 = Zpi;
% Sg1 = g1*gs1; Sg2 = g2*gs2;

%%%%
g_ypi = vrrotvec2mat([0,1,0,pi]);
% Sg1 = g_ypi*g2; Sg2 = g_ypi*g1;
% Sg1 = g_ypi*g2*Zpi; Sg2 = g_ypi*g1;
% Sg1 = g_ypi*g2;     Sg2 = g_ypi*g1*Zpi;
Sg1 = g_ypi*g2*Zpi; Sg2 = g_ypi*g1*Zpi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mvec_12_nmax = zeros(1,num_cols); SMvec_12_nmax = zeros(1,num_cols);
Mvec_ab_12  = calc_Mvec( g1,  g2, symm_orders);
SMvec_ab_12 = calc_Mvec(Sg1, Sg2, symm_orders);

mat_name = [data_fname0,'Sarr_combined_',num2str(Nmax),'.mat'];
% mat_name = [data_fname0,'Sarr_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
size(orth(S))

norm(Mvec_12_nmax*S - SMvec_12_nmax*S)

rmpath(genpath(util_dir));