clear all; clc;

addpath(genpath('../Util_functions/'));
addpath(genpath('../GB_Parameters/'));


s1 = load('rand_gb_rots.mat');
rot_mats = s1.rot_mats;

s1 = load('Sarr.mat');
Svec = s1.S;

rots = rot_mats(:,:,1);

% N = [0,4];
N = 0:4;
% N = 1;
M1 = mbp_funcs_vals(rots, N);

g_ypi = vrrotvec2mat([0,1,0, pi]);
r2 = rots(:,4:6); r1 = rots(:,1:3);
r1_ges = g_ypi*r2; r2_ges = g_ypi*r1;
M1_ges = mbp_funcs_vals([r1_ges, r2_ges], N);

norm(M1*Svec - M1_ges*Svec)

rmpath(genpath('../Util_functions/'));
rmpath(genpath('../GB_Parameters/'));