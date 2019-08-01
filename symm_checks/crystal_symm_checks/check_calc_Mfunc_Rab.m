%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Code to check extracting Mab function from Rab
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a_val = floor(rand()*6); b_val = floor(rand()*6); 
a_val = 4; b_val = 3;
nsz = (2*a_val+1)*(2*b_val+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s1 = load([top_dir,'data_files/GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;

rots1  = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
% rots1 = rot_mats(:,:,1);
rots_r = rot_mats(:,:,floor(size(rot_mats,3)*rand()));

g1 = rots1(:,1:3); g2 = rots1(:,4:6); 

ra_1 = rotmat_to_angs(g1'); check_conv(g1', ra_1);
ra_2 = rotmat_to_angs(g2'); check_conv(g2', ra_2);

R_ab_12  = so4_irrep(ra_1, ra_2, a_val, b_val);

Mvec_ab_12 = calc_Mfunc_Rab(R_ab_12, a_val, b_val);

% %%%% Check Z-rotation
% Zth = vrrotvec2mat([0,0,1,rand()*2*pi]);
% g1 = Zth*g1; g2 = Zth*g2; 
% 
% ra_1 = rotmat_to_angs(g1'); check_conv(g1', ra_1);
% ra_2 = rotmat_to_angs(g2'); check_conv(g2', ra_2);
% 
% R_ab_12  = so4_irrep(ra_1, ra_2, a_val, b_val);
% 
% Mvec1_ab_12 = calc_Mfunc_Rab(R_ab_12, a_val, b_val);
% norm(Mvec1_ab_12-Mvec_ab_12)

%%%% Add Jeremy's code Mvec
mbp_angs = rots_to_angs(g1', g2');
Mvec1_ab_12 = mbp_basis(a_val, b_val, mbp_angs);

norm(Mvec1_ab_12'-Mvec_ab_12)

rmpath(genpath(util_dir));
