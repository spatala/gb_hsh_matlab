clear all; clc;

addpath(genpath('../Util_functions/'));
addpath(genpath('../GB_Parameters/'));

s1 = load('rand_gb_rots.mat');
rot_mats = s1.rot_mats;

mat_name = 'Sarr_nmax_4.mat';
s1 = load(mat_name);
Svec = s1.Sarr_nmax;

s1 = load('O_SymmMat.mat');
SymmMat = s1.SymmMat_O;


rots = rot_mats(:,:,1);
r1 = rots(:,1:3); r2 = rots(:,4:6);

N = [0,4];
M1 = mbp_funcs_vals(rots, N);

% M1 = calc_Mfunc(a_val,b_val,rots);

% S1 = SymmMat{4}; S2 = SymmMat{12};
diff_vec = zeros(24*24,1);
ct3 = 1;
for ct1 = 1:24
    for ct2 = 1:24
        S1 = SymmMat{ct1}; S2 = SymmMat{ct2};
        
        r1s = r1*S1; r2s = r2*S2;
        M1_symm = mbp_funcs_vals([r1s,r2s], N);
        
        diff_vec(ct3) = norm(M1_symm*Svec - M1*Svec);
        ct3 = ct3 + 1;
    end
end

max(diff_vec)

rmpath(genpath('../Util_functions/'));
rmpath(genpath('../GB_Parameters/'));