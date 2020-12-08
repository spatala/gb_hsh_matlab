%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to check the right multiplication operation for MBP basis functions
%
% Mvec_ab_12 = M(g1,g2)
% gt1 = g1*(gr1'); gt2 = g2*(gr2');
% tMvec_ab_12 = M(gt1,gt2)
% rMat1 = (kron(U_a (gr1), U_b (gr2)).'); 
% rMat = kron(eye(ng), rMat1);
% Check that tMvec_ab_12 == Mvec_ab_12*(rMat))
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
top_dir = get_top_dir();
util_dir = strcat(top_dir,'Util_functions','/'); 
addpath(genpath(util_dir));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_val = floor(rand()*6); b_val = floor(rand()*6); 
nsz = (2*a_val+1)*(2*b_val+1);
c_val = min(a_val, b_val); ng = (2*c_val + 1); nsz1 = nsz*ng;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s1 = load([top_dir,'data_files/GB_Parameters/rand_gb_rots.mat']); 
rot_mats = s1.rot_mats;
num = 1000;

diff_vec = zeros(num,1);

for tct1=1:num
rots1  = rot_mats(:,:,1+floor(size(rot_mats,3)*rand()));
rots_r = rot_mats(:,:,1+floor(size(rot_mats,3)*rand()));

Zpi = vrrotvec2mat([0,0,1,rand()*2*pi]);
g1 = Zpi*rots1(:,1:3); g2 = Zpi*rots1(:,4:6); 
Zpi = vrrotvec2mat([0,0,1,rand()*2*pi]);
gr1 = Zpi*rots_r(:,1:3); gr2 = Zpi*rots_r(:,4:6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get Mvecs
mbp_angs = rots_to_angs(g1, g2);
% .' - for simple transpose and not conjugate transpose
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

rMat = kron(eye(ng), rMat1);
if (max(abs(tMvec_ab_12 - Mvec_ab_12*(rMat))) > 1e-8)
    disp(tct1)
    disp([tct1,[ra_r1, ra_r2]*180/pi])
    disp(max(abs(tMvec_ab_12 - Mvec_ab_12*(rMat))))
end
diff_vec(tct1) = max(abs(tMvec_ab_12 - Mvec_ab_12*(rMat)));

end
display(max(diff_vec))

rmpath(genpath(util_dir));