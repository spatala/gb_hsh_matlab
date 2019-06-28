%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Code to check the rotation of SO(4) Functions R^(a,b)
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
addpath(genpath(top_dir));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_val = floor(rand()*10); Na = 2*a_val;
b_val = floor(rand()*10); Nb = 2*b_val;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;

rots1  = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
rots_r = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
rots_l = rot_mats(:,:,floor(size(rot_mats,3)*rand())); 

g1  = rots1(:,1:3);  g2  = rots1(:,4:6); 
gr1 = rots_r(:,1:3); gr2 = rots_r(:,4:6);
gl1 = rots_l(:,1:3); gl2 = rots_l(:,4:6); 

R_ab_12  = so4_irrep(g1 ,g2 ,Na,Nb);
Rr_ab_12 = so4_irrep(gr1,gr2,Na,Nb);
Rl_ab_12 = so4_irrep(gl1,gl2,Na,Nb);


rot1 = Rr_ab_12*R_ab_12*Rl_ab_12;
tmat1 = gl1*g1*gr1; tmat2 = gl2*g2*gr2;
R_ab_t1t2 = so4_irrep(tmat1,tmat2,Na,Nb);
norm(R_ab_t1t2 - rot1)

%%%%% Check the tranpose rotations as well!
trR_ab_12  = transpose(R_ab_12);
trRr_ab_12 = transpose(Rr_ab_12);
trRl_ab_12 = transpose(Rl_ab_12);
tr_rot1 = trRl_ab_12*trR_ab_12*trRr_ab_12;
trR_ab_t1t2 = transpose(R_ab_t1t2);
norm(trR_ab_t1t2 - tr_rot1)
