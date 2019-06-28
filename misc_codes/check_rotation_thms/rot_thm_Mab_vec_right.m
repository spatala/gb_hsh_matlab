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
% a_val = floor(rand()*6); b_val = floor(rand()*6); 
a_val = 4; b_val = 3;
c_val = min(a_val, b_val);
Na = 2*a_val; Nb = 2*b_val; Nc = 2*c_val; nsz = (Na+1)*(Nb+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;

% rots1  = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
rots1 = rot_mats(:,:,1);
rots_r = rot_mats(:,:,floor(size(rot_mats,3)*rand()));


g1 = rots1(:,1:3); g2 = rots1(:,4:6); 
% gr1 = rots_r(:,1:3); gr2 = rots_r(:,4:6); 
gr1 = vrrotvec2mat([0,0,1,pi]); gr2 = eye(3); 
% gr1 = eye(3,3); gr2 = eye(3,3);

R_ab_12  = so4_irrep(g1 ,g2 ,Na,Nb);
Rr_ab_12 = so4_irrep(gr1,gr2,Na,Nb);

trR_ab_12 = transpose(R_ab_12);
Mvec_ab_12 = calc_Mfunc_Rab(R_ab_12, a_val, b_val);

trRr_ab_12 = transpose(Rr_ab_12); 

Rmult = Rr_ab_12*R_ab_12;
Mvec1_ab_12 = calc_Mfunc_Rab(Rmult, a_val, b_val);

tr_Rmult = transpose(Rmult);
tr_Rmult_v1 = transpose(Rmult(:));

Rvec_ab_12 = reshape(transpose(R_ab_12), [1,nsz*nsz]);
trRv_ab_12 = transpose(R_ab_12(:));


tr_Rrot = kron(eye(size(Rr_ab_12)),trRr_ab_12);
tr_Rmult_v2 = trRv_ab_12*tr_Rrot;

Mrot = calc_Mrot_mat_Rab(tr_Rrot, a_val, b_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mrot1 = kron(eye(Nc+1), trRr_ab_12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% norm(Mrot - Mrot1)

Mvec2_ab_12 = Mvec_ab_12*Mrot1;
% norm(tr_Rmult_v1 - tr_Rmult_v2)
norm(Mvec1_ab_12-Mvec2_ab_12)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mvec1_ab_12 = calc_Mfunc_Rab(Rmult, a_val, b_val);

% tR_ab_12 = so4_irrep(gl1*g1*gr1,gl2*g2*gr2,Na,Nb);
% Mvec2_ab_12 = calc_Mfunc_Rab_full(Rmult, a_val, b_val);
% norm(Mvec1_ab_12 - Mvec2_ab_12)

% Mmult_v1 = calc_Mfunc_Rab(Rmult, a_val, b_val);
% Rrot = kron(transpose(Rr_ab_12),Rl_ab_12);
% tr_Rrot = transpose(Rrot);
% Mrot = calc_Mrot_mat_Rab(Rrot, a_val, b_val);
% norm(Mmult_v1 - Mvec_ab_12*Mrot)
% Mrot(abs(Mrot)<1e-10) = 0;
% Mrot = sparse(Mrot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmpath(genpath(util_dir));