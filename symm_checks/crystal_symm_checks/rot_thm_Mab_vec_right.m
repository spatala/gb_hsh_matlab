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
nsz = (2*a_val+1)*(2*b_val+1);
c_val = min(a_val, b_val); ng = (2*c_val + 1); nsz1 = nsz*ng;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s1 = load([top_dir,'data_files/GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;
num = 1000;

diff_vec = zeros(num,1);

for tct1=1:num
rots1  = rot_mats(:,:,1+floor(size(rot_mats,3)*rand()));
% rots1 = rot_mats(:,:,1);
rots_r = rot_mats(:,:,1+floor(size(rot_mats,3)*rand()));

Zpi = vrrotvec2mat([0,0,1,rand()*2*pi]);
g1 = Zpi*rots1(:,1:3); g2 = Zpi*rots1(:,4:6); 
Zpi = vrrotvec2mat([0,0,1,rand()*2*pi]);
gr1 = Zpi*rots_r(:,1:3); gr2 = Zpi*rots_r(:,4:6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get Mvecs
mbp_angs = rots_to_angs(g1, g2);
Mvec_ab_12 = (mbp_basis(a_val, b_val, mbp_angs).');

gt1 = g1*(gr1'); gt2 = g2*(gr2');
tmbp_angs = rots_to_angs(gt1, gt2);
tMvec_ab_12 = (mbp_basis(a_val, b_val, tmbp_angs).');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ra_r1 = rotmat_to_angs(gr1); check_conv(gr1, ra_r1);
ra_r2 = rotmat_to_angs(gr2); check_conv(gr2, ra_r2);

U_a = rotation_mat(a_val, ra_r1);
U_b = rotation_mat(b_val, ra_r2);
rMat1 = (kron(U_a, U_b).');

% %%%% Check on the nsz components of Mvec
% tn1 = floor(rand()*ng); inds = (tn1*nsz + 1):((tn1+1)*nsz);
% Mvec1 = Mvec_ab_12(inds); tMvec1 = tMvec_ab_12(inds); 
% if (norm(tMvec1 - Mvec1*(rMat1)) > 1e-10)
%     disp(tct1)
% else
%     disp(max(abs(tMvec1 - Mvec1*(rMat1))))
% end

rMat = kron(eye(ng), rMat1);
if (max(abs(tMvec_ab_12 - Mvec_ab_12*(rMat))) > 1e-10)
    disp(tct1)
    disp([tct1, max(abs(tMvec_ab_12 - Mvec_ab_12*(rMat))),[ra_r1, ra_r2]*180/pi])
end
diff_vec(tct1) = max(abs(tMvec_ab_12 - Mvec_ab_12*(rMat)));

end

rmpath(genpath(util_dir));