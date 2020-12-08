%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to check the rotation operations of SO(3) Functions U^a
% $U(gl)*U(g)*U(gr) = U(gl*g1*gr)$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
top_dir = get_top_dir();
util_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(top_dir));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = load([top_dir,'data_files/GB_Parameters/rand_gb_rots.mat']);
rot_mats = s1.rot_mats;
rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
r1 = rots1(:,1:3);
a_val = floor(rand()*12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rot_angs = rotmat_to_angs(r1); check_conv(r1, rot_angs);
U_a = rotation_mat(a_val, rot_angs);
trots = rot_mats(:,:,floor(size(rot_mats,3)*rand())); 
tl = trots(:,1:3); tr = trots(:,4:6);
ra_l = rotmat_to_angs(tl); check_conv(tl, ra_l);
ra_r = rotmat_to_angs(tr); check_conv(tr, ra_r);
U_l = rotation_mat(a_val, ra_l);
U_r = rotation_mat(a_val, ra_r);
tmat = tl*r1*tr;
ra_t = rotmat_to_angs(tmat); check_conv(tmat, ra_t);
U_t = rotation_mat(a_val, ra_t);
norm(U_l*U_a*U_r - U_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rmpath(genpath(util_dir));
