clear all; clc;

curr_pwd = split(pwd,'/');
top_dir = '';
for ct1=1:length(curr_pwd)
    top_dir = strcat(top_dir,curr_pwd{ct1},'/');
    if (strcmp(curr_pwd{ct1},'gb_hsh_matlab'))
        break;
    end
end
util_dir = strcat(top_dir,'Util_functions/');
gb_data_dir = strcat(top_dir,'GB_Parameters/');
addpath(genpath(util_dir)); addpath(genpath(gb_data_dir));

s1 = set_vars();
Nmax = s1.Nmax; pt_grp = s1.pt_grp;
fname = [top_dir,'data_files', '/ptgrp_',pt_grp];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mat_name = [fname,'/cryst_symm/symm_ab_',...
pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
symm_orders = s1.symm_orders;
nsymm = size(symm_orders,1);
a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1, b1);
num_rows = sum((2*a1+1).*(2*b1+1).*(2*c1+1));


s1 = load('rand_gb_rots.mat');
rot_mats = s1.rot_mats;
rots = rot_mats(:,:,1);
r1 = rots(:,1:3); r2 = rots(:,4:6);

M1 = calc_Mfunc_ab_array(symm_orders, num_rows, rots);

mat_name = [fname,'/ge_symm/Sarr_ges_Nmax_',...
    num2str(Nmax),'.mat'];
s1 = load(mat_name);
Svec = s1.S;

% %%%% Doesn't work!! (as expected)
% mat_name = [fname,'/ge_symm/orthY_ges_Nmax_',...
%     num2str(Nmax),'.mat'];
% s1 = load(mat_name); Svec = s1.Y1;
% mat_name = [fname,'/ge_symm/Y_ges_Nmax_',...
%     num2str(Nmax),'.mat'];
% s1 = load(mat_name); Svec = s1.col1;


g_ypi = vrrotvec2mat([0,1,0, pi]);
r1_ges = g_ypi*r2; r2_ges = g_ypi*r1;
rots_ges = [r1_ges, r2_ges];
M1_ges = calc_Mfunc_ab_array(symm_orders, num_rows, rots_ges);

norm(M1*Svec - M1_ges*Svec)

rmpath(genpath(util_dir)); rmpath(genpath(gb_data_dir));