clear all; clc;

curr_pwd = split(pwd,'/');
top_dir = '';
for ct1=1:length(curr_pwd)
    top_dir = strcat(top_dir,curr_pwd{ct1},'/');
    if (strcmp(curr_pwd{ct1},'MATLAB_Codes'))
        break;
    end
end
util_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(top_dir));

s1 = set_vars();
Nmax = s1.Nmax; pt_grp = s1.pt_grp;

data_fname0 = [top_dir,'data_files/ptgrp_',pt_grp,'/cryst_symm/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fname1 = [top_dir,'data_files/ptgrp_',pt_grp];
symm_mat_name = [fname1,'/SymmMat_',pt_grp,'.mat'];
s1 = load(symm_mat_name);
SymmMat = s1.SymmMat;
num_symm = length(SymmMat);

mat_name = [fname1,'/cryst_symm/symm_ab_',...
    pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
symm_orders = s1.symm_orders;
nsymm = size(symm_orders,1);
a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1,b1);
num_rows = sum((2*a1+1).*(2*b1+1).*(2*c1+1));

s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;
rots = rot_mats(:,:,1); r1 = rots(:,1:3); r2 = rots(:,4:6);

mat_name = [fname1,'/cryst_symm/Sarr_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
Svec = s1.S;

M1 = calc_Mfunc_ab_array(symm_orders, num_rows, rots);

diff_vec = zeros(num_symm*num_symm,1);
ct3 = 1;
for ct1 = 1:num_symm
    for ct2 = 1:num_symm
        S1 = SymmMat{ct1}; S2 = SymmMat{ct2};
        
        r1s = r1*S1; r2s = r2*S2;
        rots_symm = [r1s, r2s];
        M1_symm = calc_Mfunc_ab_array(symm_orders, num_rows, rots_symm);
        
        diff_vec(ct3) = norm(M1_symm*Svec - M1*Svec);
        ct3 = ct3 + 1;
    end
end

max(diff_vec)

rmpath(genpath(top_dir));