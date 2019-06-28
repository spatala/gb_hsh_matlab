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
Na = 2*a_val; Nb = 2*b_val; nsz = (Na+1)*(Nb+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;

rots1  = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
rots_l = rot_mats(:,:,floor(size(rot_mats,3)*rand())); 

g1  = rots1(:,1:3);  g2  = rots1(:,4:6); 
gl1 = rots_l(:,1:3); gl2 = rots_l(:,4:6); 

R_ab_12  = so4_irrep(g1 ,g2 ,Na,Nb);
Rv_ab_12 = reshape(transpose(R_ab_12),[1,nsz*nsz]);

Rl_ab_12 = so4_irrep(gl1,gl2,Na,Nb);
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
