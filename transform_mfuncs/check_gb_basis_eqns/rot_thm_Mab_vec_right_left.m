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
a_val = floor(rand()*6); b_val = floor(rand()*6); c_val = min(a_val, b_val);
% a_val = 3; b_val = 3; c_val = min(a_val, b_val);
Na = 2*a_val; Nb = 2*b_val; Nc = 2*c_val; nsz = (Na+1)*(Nb+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;

rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand())); 
g1 = rots1(:,1:3); g2 = rots1(:,4:6); 
ax_ang_1 = vrrotmat2vec(g1); ax_ang_2 = vrrotmat2vec(g2);
U1  = rotation( ax_ang_1(1:3),  ax_ang_1(4), Na);
U2  = rotation( ax_ang_2(1:3),  ax_ang_2(4), Nb);

rots_r = rot_mats(:,:,floor(size(rot_mats,3)*rand())); 
gr1 = rots_r(:,1:3); gr2 = rots_r(:,4:6); 
% gr1 = eye(3,3); gr2 = eye(3,3);
ax_ang_r1 = vrrotmat2vec(gr1); ax_ang_r2 = vrrotmat2vec(gr2);
Ur1 = rotation(ax_ang_r1(1:3), ax_ang_r1(4), Na);
Ur2 = rotation(ax_ang_r2(1:3), ax_ang_r2(4), Nb);

% rots_l = rot_mats(:,:,floor(size(rot_mats,3)*rand())); 
% gl1 = rots_l(:,1:3); gl2 = rots_l(:,4:6); 
% gl1 = eye(3,3); gl2 = eye(3,3);
%%%%% Only works (so far) for 
%%%%% 1) X,Y,Z-pi rotations!; 2) Z-theta rotations!
Ypi = [0,1,0,pi]; gl1 = vrrotvec2mat(Ypi); gl2 = vrrotvec2mat(Ypi);
% Zth = [0,0,1,rand()]; gl1 = vrrotvec2mat(Zth); gl2 = vrrotvec2mat(Zth);
ax_ang_l1 = vrrotmat2vec(gl1); ax_ang_l2 = vrrotmat2vec(gl2);
Ul1 = rotation(ax_ang_l1(1:3), ax_ang_l1(4), Na);
Ul2 = rotation(ax_ang_l2(1:3), ax_ang_l2(4), Nb);

R_ab_12 = kron(U1, U2); Rr_ab_12 = kron(Ur1, Ur2); Rl_ab_12 = kron(Ul1, Ul2);
Mvec_ab_12 = calc_Mfunc_Rab(R_ab_12, a_val, b_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rmult = Rl_ab_12*R_ab_12*Rr_ab_12;
Mmult_v1 = calc_Mfunc_Rab(Rmult, a_val, b_val);
Rrot = kron(transpose(Rl_ab_12),Rr_ab_12);
Mrot = calc_rot_mat_Rab(Rrot, a_val, b_val);
norm(Mmult_v1 - Mvec_ab_12*Mrot)
% Mrot(abs(Mrot)<1e-10) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%