clear all; clc;
clear all; clc;


fname = get_dir_name();
Nmax = 4;
% pt_grp = 'O';
pt_grp = 'C2';

mat_name = [fname,'/ptgrp_',pt_grp,'/cryst_symm/symm_ab_',...
    pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
symm_orders = s1.symm_orders;
nsymm = size(symm_orders,1);
a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1, b1);
num_rows = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
Nmax = max(a1);

s1 = load('rand_gb_rots.mat');
rot_mats = s1.rot_mats;

mat_name = [fname,'/ptgrp_',pt_grp,'/ge_symm/Sarr_ges_Nmax_',...
    num2str(Nmax),'.mat'];
s1 = load(mat_name);
Svec = s1.S;
% s1 = load('Y_ges.mat');
% Svec = orth(full(s1.col1));

rots = rot_mats(:,:,1);
r1 = rots(:,1:3); r2 = rots(:,4:6);

M1 = calc_Mfunc_ab_array(symm_orders, num_rows, rots);

g_ypi = vrrotvec2mat([0,1,0, pi]);
r1_ges = g_ypi*r2; r2_ges = g_ypi*r1;
rots_ges = [r1_ges, r2_ges];
M1_ges = calc_Mfunc_ab_array(symm_orders, num_rows, rots_ges);

norm(M1*Svec - M1_ges*Svec)
