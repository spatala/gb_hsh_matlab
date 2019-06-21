clear all; clc;

addpath(genpath('../Util_functions/'));
addpath(genpath('../GB_Parameters/'));

fname = get_dir_name();

% pt_grp = 'O';
pt_grp = 'C2';
Nmax = 4;
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


rots = rot_mats(:,:,1);
r1 = rots(:,1:3); r2 = rots(:,4:6);



M1 = calc_Mfunc_ab_array(symm_orders, num_rows, rots);

% a_val = 4; b_val = 4;
% M1 = mbp_funcs_vals(rots, N);
% M1 = calc_Mfunc(a_val,b_val,rots);
% norm(M1*Svec - M1_ges*Svec)

mat_name = [fname,'/ptgrp_',pt_grp,'/combined_symm/Sarr_combined_Nmax_',...
    num2str(Nmax),'.mat'];
s1 = load(mat_name);
Svec = s1.S;

symm_mat_name = [fname,'/ptgrp_',pt_grp,'/SymmMat_',pt_grp,'.mat'];
s1 = load(symm_mat_name);
SymmMat = s1.SymmMat;
num_symm = length(SymmMat);

diff_vec = zeros(num_symm*num_symm*2,1);
ct3 = 1;
for ct1 = 1:num_symm
    for ct2 = 1:num_symm
        S1 = SymmMat{ct1}; S2 = SymmMat{ct2};

        r1s = r1*S1; r2s = r2*S2;
        rots_symm = [r1s, r2s];
        M1_symm = calc_Mfunc_ab_array(symm_orders, num_rows, rots_symm);

        diff_vec(ct3) = norm(M1_symm*Svec - M1*Svec);
        ct3 = ct3 + 1;
        
        g_ypi = vrrotvec2mat([0,1,0, pi]);
        r1_ges = g_ypi*r2s; r2_ges = g_ypi*r1s;
        rots_ges = [r1_ges, r2_ges];
        M1_ges = calc_Mfunc_ab_array(symm_orders, num_rows, rots_ges);
        diff_vec(ct3) = norm(M1_ges*Svec - M1*Svec);
        ct3 = ct3 + 1;
    end
end
max(diff_vec)

rmpath(genpath('../Util_functions/'));
rmpath(genpath('../GB_Parameters/'));