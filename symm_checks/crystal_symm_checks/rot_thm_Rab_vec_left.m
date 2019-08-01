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
addpath(genpath(util_dir));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_val = floor(rand()*6); b_val = floor(rand()*6); 
% a_val = 3; b_val = 2;
nsz = (2*a_val+1)*(2*b_val+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s1 = load([top_dir,'data_files/GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;

rots1  = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
rots_l = rot_mats(:,:,floor(size(rot_mats,3)*rand())); 

g1  = rots1(:,1:3);  g2  = rots1(:,4:6);
ra_1 = rotmat_to_angs(g1'); check_conv(g1', ra_1);
ra_2 = rotmat_to_angs(g2'); check_conv(g2', ra_2);

gl1 = rots_l(:,1:3); gl2 = rots_l(:,4:6); 
ra_l1 = rotmat_to_angs(gl1'); check_conv(gl1', ra_l1);
ra_l2 = rotmat_to_angs(gl2'); check_conv(gl2', ra_l2);


R_ab_12  = so4_irrep(ra_1, ra_2, a_val, b_val);
Rv_ab_12 = reshape(transpose(R_ab_12), [1,nsz*nsz]);

Rl_ab_12 = so4_irrep(ra_l1, ra_l2, a_val, b_val);

Rmult = R_ab_12*Rl_ab_12;
Rmult_v1 = reshape(transpose(Rmult),[1,nsz*nsz]);
Rmult_v2 = Rv_ab_12*kron(eye(nsz,nsz),Rl_ab_12);

norm(Rmult_v1 - Rmult_v2)

%%%%% Check the tranpose rotations as well!
trR_ab_12 = transpose(R_ab_12);
trRv_ab_12 = transpose(R_ab_12(:));
trRl_ab_12 = transpose(Rl_ab_12);
tr_Rmult = trRl_ab_12*trR_ab_12;
tr_Rmult_v1 = transpose(Rmult(:));
tr_Rmult_v2 = trRv_ab_12*kron(Rl_ab_12,eye(nsz,nsz));
norm(tr_Rmult_v1 - tr_Rmult_v2)


rmpath(genpath(util_dir));
