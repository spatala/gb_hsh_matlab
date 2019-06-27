%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Code to check the rotation of SO(4) Functions R^(a,b)
%%%%%%%
clear all; clc;

% check_symm='ypi';
% check_symm='flip';
check_symm='ges';

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
% Nmax = 1;
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
if (strcmp(check_symm, 'ypi'))
    gypi = vrrotvec2mat([0,1,0,pi]); g1_s = gypi*g1; g2_s = gypi*g2;
    symm_mat = generate_ypi_left(symm_orders);
end
if (strcmp(check_symm, 'flip'))
    g1_s = g2; g2_s = g1;
    symm_mat = generate_flip_mat(symm_orders);
end
if (strcmp(check_symm, 'ges'))
    gypi = vrrotvec2mat([0,1,0,pi]); 
    g1_s = gypi*g2; g2_s = gypi*g1;

    Ypi_symm_mat = generate_ypi_left(symm_orders);
    flip_mat = generate_flip_mat(symm_orders);
    symm_mat = Ypi_symm_mat*flip_mat;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mvec_12_nmax = zeros(1,num_cols);
SMvec_12_nmax = zeros(1,num_cols);
ind_start = 1;

for ct1=1:nsymm
    a_val = symm_orders(ct1,1); b_val = symm_orders(ct1,2);
    c_val = min(a_val, b_val);
    Na = 2*a_val; Nb = 2*b_val; Nc = 2*c_val; nsz = (Na+1)*(Nb+1);
    
    R_ab_12 = so4_irrep(g1,g2,Na,Nb);
    Mvec_ab_12 = calc_Mfunc_Rab(R_ab_12, a_val, b_val);
    
    Ry_ab_12 = so4_irrep(g1_s,g2_s,Na,Nb);
    SMvec_ab_12 = calc_Mfunc_Rab(Ry_ab_12, a_val, b_val);
    
    nsz1 = nsz*(Nc+1);
    ind_stop = ind_start + nsz1 - 1;
    % [ind_start, ind_stop]
    Mvec_12_nmax(ind_start:ind_stop) = Mvec_ab_12;
    SMvec_12_nmax(ind_start:ind_stop) = SMvec_ab_12;
    ind_start = ind_stop + 1;
end

norm(Mvec_12_nmax*symm_mat - SMvec_12_nmax)
[v, d] = eig(full(symm_mat));
col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
S = orth(v(:,col));
mat_name = [data_fname0,'ges_nmax_',num2str(Nmax),'.mat'];
save(mat_name, 'S');

% s1 = load(mat_name); S = s1.S;
% size(orth(S))

norm(Mvec_12_nmax*S - SMvec_12_nmax*S)

rmpath(genpath(util_dir));