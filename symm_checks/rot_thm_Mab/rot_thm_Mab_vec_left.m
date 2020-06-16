%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to check the left multiplication operation (by Ypi) 
% for MBP basis functions
%
% Mvec_ab_12 = M(g1,g2)
% gt1 = Ypi*g1; gt2 = Ypi*g2;
% tMvec_ab_12 = M(gt1,gt2)
% rMat = generate_ypi_left_ab(a_val, b_val)
% Check that tMvec_ab_12 == Mvec_ab_12*(rMat))
%

clear all; clc;
top_dir = get_top_dir();
util_dir = strcat(top_dir,'Util_functions','/'); 
addpath(genpath(util_dir));
symmgen_dir = strcat(top_dir,'symm_checks','/'); 
addpath(genpath(symmgen_dir));

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

Zth = vrrotvec2mat([0,0,1,rand()*2*pi]);
g1 = Zth*rots1(:,1:3); g2 = Zth*rots1(:,4:6); 
Ypi = vrrotvec2mat([0,1,0,pi]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get Mvecs
mbp_angs = rots_to_angs(g1, g2);
Mvec_ab_12 = (mbp_basis(a_val, b_val, mbp_angs)');

gt1 = Ypi*g1; gt2 = Ypi*g2;
tmbp_angs = rots_to_angs(gt1, gt2);
tMvec_ab_12 = (mbp_basis(a_val, b_val, tmbp_angs)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rMat = generate_ypi_left_ab(a_val, b_val);

diff_vec(tct1) = max(abs(tMvec_ab_12 - Mvec_ab_12*(rMat)));

end

display(max(diff_vec))

rmpath(genpath(util_dir));
rmpath(genpath(symmgen_dir));