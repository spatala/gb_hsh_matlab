%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to check the multiplication operations of SO(4) Functions R^(a,b)
% Rr_ab_12 = R(gr1',gr2')
% Rl_ab_12 = eye(nsz,nsz)
% R_ab_12 = R(g1',g2')
% Rv_ab_12 = reshape(transpose(R_ab_12), [1,nsz*nsz])
% 
% Rmult = Rr_ab_12*R_ab_12*Rl_ab_12;
% Rmult_v1 = reshape(transpose(Rmult),[1,nsz*nsz]);
% Rmult_v2 = Rv_ab_12*kron(transpose(Rr_ab_12),Rl_ab_12);
% 
% Show that Rmult_v1 == Rmult_v2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;
top_dir = get_top_dir();
util_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(top_dir));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_val = floor(rand()*6); b_val = floor(rand()*6); 
% a_val = 3; b_val = 2;
nsz = (2*a_val+1)*(2*b_val+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s1 = load([top_dir,'data_files/GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;
rots1  = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
rots_r = rot_mats(:,:,floor(size(rot_mats,3)*rand()));

g1  = rots1(:,1:3);  g2  = rots1(:,4:6);
ra_1 = rotmat_to_angs(g1'); check_conv(g1', ra_1);
ra_2 = rotmat_to_angs(g2'); check_conv(g2', ra_2);

gr1 = rots_r(:,1:3); gr2 = rots_r(:,4:6);
ra_r1 = rotmat_to_angs(gr1'); check_conv(gr1', ra_r1);
ra_r2 = rotmat_to_angs(gr2'); check_conv(gr2', ra_r2);


R_ab_12  = so4_irrep(ra_1, ra_2, a_val, b_val);
Rv_ab_12 = reshape(transpose(R_ab_12), [1,nsz*nsz]);

Rr_ab_12 = so4_irrep(ra_r1, ra_r2, a_val, b_val);

Rmult = Rr_ab_12*R_ab_12;
Rmult_v1 = reshape(transpose(Rmult),[1,nsz*nsz]);
Rmult_v2 = Rv_ab_12*kron(transpose(Rr_ab_12),eye(nsz,nsz));

norm(Rmult_v1 - Rmult_v2)

%%%%% Check the tranpose rotations as well!
trR_ab_12 = transpose(R_ab_12);
trRv_ab_12 = transpose(R_ab_12(:));
trRr_ab_12 = transpose(Rr_ab_12);
tr_Rmult = trR_ab_12*trRr_ab_12;
tr_Rmult_v1 = transpose(Rmult(:));
tr_Rmult_v2 = trRv_ab_12*kron(eye(nsz,nsz),trRr_ab_12);
norm(tr_Rmult_v1 - tr_Rmult_v2)

rmpath(genpath(util_dir));