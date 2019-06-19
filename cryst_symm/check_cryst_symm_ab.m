clear all; clc;

addpath(genpath('../Util_functions/'));
addpath(genpath('../GB_Parameters/'));

s1 = load('rand_gb_rots.mat'); rot_mats = s1.rot_mats;

% pt_grp = 'O';
pt_grp = 'C2';

Nmax = 4;
fname = get_dir_name();
symm_mat_name = [fname,'/ptgrp_',pt_grp,'/SymmMat_',pt_grp,'.mat'];
s1 = load(symm_mat_name);
SymmMat = s1.SymmMat;
num_symm = length(SymmMat);

rots = rot_mats(:,:,1);
r1 = rots(:,1:3); r2 = rots(:,4:6);
for a_val = 0:Nmax
    for b_val = 0:Nmax
        
        mat_name = [fname,'/ptgrp_',pt_grp,'/cryst_symm/Sarr_',num2str(a_val),'_',num2str(b_val),'.mat'];
        s1 = load(mat_name);
        Svec = s1.S;
        
        
        M1 = calc_Mfunc(a_val,b_val,rots);
        
        diff_vec = zeros(num_symm*num_symm,1);
        ct3 = 1;
        for ct1 = 1:num_symm
            for ct2 = 1:num_symm
                S1 = SymmMat{ct1}; S2 = SymmMat{ct2};
                
                r1s = r1*S1; r2s = r2*S2;
                M1_symm = calc_Mfunc(a_val,b_val,[r1s, r2s]);
                
                diff_vec(ct3) = norm(M1_symm*Svec - M1*Svec);
                ct3 = ct3 + 1;
            end
        end
        [a_val, b_val, size(Svec,2),max(diff_vec)]
    end
end
rmpath(genpath('../Util_functions/'));
rmpath(genpath('../GB_Parameters/'));