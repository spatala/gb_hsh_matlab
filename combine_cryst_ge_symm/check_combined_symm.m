clear all; clc;

addpath(genpath('./Util_functions/'));
addpath(genpath('./GB_Parameters/'));


s1 = load('rand_gb_rots.mat');
rot_mats = s1.rot_mats;

s1 = load('Sarr_C2_4_4_combined.mat');
Svec = s1.S;

rots = rot_mats(:,:,1);
r1 = rots(:,1:3); r2 = rots(:,4:6);

a_val = 4; b_val = 4;
% M1 = mbp_funcs_vals(rots, N);
M1 = calc_Mfunc(a_val,b_val,rots);


% norm(M1*Svec - M1_ges*Svec)


fname = get_dir_name();
pt_grp = 'C2';
symm_mat_name = [fname,'/ptgrp_',pt_grp,'/SymmMat_',pt_grp,'.mat'];
s1 = load(symm_mat_name);
SymmMat = s1.SymmMat;
num_symm = length(SymmMat);

% s1 = load('data_files/ptgrp_C2/cryst_symm/Sarr_4_4.mat');
% Svec = s1.S;


diff_vec = zeros(num_symm*num_symm*2,1);
ct3 = 1;
for ct1 = 1:num_symm
    for ct2 = 1:num_symm
        S1 = SymmMat{ct1}; S2 = SymmMat{ct2};

        r1s = r1*S1; r2s = r2*S2;
        M1_symm = calc_Mfunc(a_val,b_val,[r1s, r2s]);

        diff_vec(ct3) = norm(M1_symm*Svec - M1*Svec);
        ct3 = ct3 + 1;
        
        g_ypi = vrrotvec2mat([0,1,0, pi]);
        r1_ges = g_ypi*r2s; r2_ges = g_ypi*r1s;
        M1_ges = calc_Mfunc(a_val,b_val,[r1_ges, r2_ges]);
        diff_vec(ct3) = norm(M1_symm*Svec - M1*Svec);
        ct3 = ct3 + 1;
    end
end
max(diff_vec)

rmpath(genpath('../Util_functions/'));
rmpath(genpath('../GB_Parameters/'));