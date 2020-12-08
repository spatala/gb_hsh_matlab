%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to check that the calculation of MBP basis function for a fixed 
% (a,b), M^{a,b}, with rotations (g1, g2) is exactly the same as the basis
% function with rotation (Z*g1, Z*g2), where Z is a random rotation along
% the z-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
top_dir = get_top_dir();
util_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(util_dir));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_val = floor(rand()*10); b_val = floor(rand()*10); 
% display([a_val, b_val])
% a_val = 4; b_val = 3;
nsz = (2*a_val+1)*(2*b_val+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s1 = load([top_dir,'data_files/GB_Parameters/rand_gb_rots.mat']); 
rot_mats = s1.rot_mats;

rots1  = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
rots_r = rot_mats(:,:,floor(size(rot_mats,3)*rand()));

g1 = rots1(:,1:3); g2 = rots1(:,4:6); 

mbp_angs = rots_to_angs(g1, g2);
Mvec_ab_12 = mbp_basis(a_val, b_val, mbp_angs);

%%%% Check Z-rotation
Zth = vrrotvec2mat([0,0,1,rand()*2*pi]);
g1z = Zth*g1; g2z = Zth*g2; 

mbp_angs_z = rots_to_angs(g1z, g2z);
Mvec1_ab_12 = mbp_basis(a_val, b_val, mbp_angs_z);

norm(Mvec1_ab_12-Mvec_ab_12)

rmpath(genpath(util_dir));
