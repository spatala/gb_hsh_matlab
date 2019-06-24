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
s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;
rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand())); 
g1 = rots1(:,1:3); ax_ang_1 = vrrotmat2vec(g1);
g2 = rots1(:,4:6); ax_ang_2 = vrrotmat2vec(g2);
rots_r = rot_mats(:,:,floor(size(rot_mats,3)*rand())); 
gr1 = rots_r(:,1:3); ax_ang_r1 = vrrotmat2vec(gr1);
gr2 = rots_r(:,4:6); ax_ang_r2 = vrrotmat2vec(gr2);
rots_l = rot_mats(:,:,floor(size(rot_mats,3)*rand())); 
gl1 = rots_l(:,1:3); ax_ang_l1 = vrrotmat2vec(gl1);
gl2 = rots_l(:,4:6); ax_ang_l2 = vrrotmat2vec(gl2);


a_val = floor(rand()*10); Na = 2*a_val;
b_val = floor(rand()*10); Nb = 2*b_val;

U1  = rotation( ax_ang_1(1:3),  ax_ang_1(4), Na);
Ur1 = rotation(ax_ang_r1(1:3), ax_ang_r1(4), Na);
Ul1 = rotation(ax_ang_l1(1:3), ax_ang_l1(4), Na);

U2  = rotation( ax_ang_2(1:3),  ax_ang_2(4), Nb);
Ur2 = rotation(ax_ang_r2(1:3), ax_ang_r2(4), Nb);
Ul2 = rotation(ax_ang_l2(1:3), ax_ang_l2(4), Nb);

R_ab_12 = kron(U1, U2);
Rr_ab_12 = kron(Ur1, Ur2);
Rl_ab_12 = kron(Ul1, Ul2);

rot1 = Rl_ab_12*R_ab_12*Rr_ab_12;

tmat1 = gl1*g1*gr1; tax_ang_1 = vrrotmat2vec(tmat1);
tmat2 = gl2*g2*gr2; tax_ang_2 = vrrotmat2vec(tmat2);

tU1 = rotation(tax_ang_1(1:3), tax_ang_1(4), Na);
tU2 = rotation(tax_ang_2(1:3), tax_ang_2(4), Nb);
R_ab_t1t2 = kron(tU1,tU2);

norm(R_ab_t1t2 - rot1)
