%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code to create random set of coefficients and check the GES
%%%%% relationship (Not complete yet!).

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

s1 = set_vars();
Nmax = s1.Nmax; pt_grp = s1.pt_grp;

data_fname0 = [top_dir,'data_files/ptgrp_',pt_grp,'/cryst_symm/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mat_name = ['Y_ges_Nmax_',num2str(Nmax),'.mat'];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s1 = load(mat_name); col1 = s1.col1;
% 
% nvec = size(col1,2);
% coeffs = col1*(rand(nvec,1)*100);




% symm_orders = zeros(Nmax^2,2);
% ct3 = 1;
% for ct1=0:Nmax
%     for ct2=0:Nmax
%         symm_orders(ct3,:) = [ct1,ct2];
%         ct3 = ct3 + 1;
%     end
% end
% 
% a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1,b1);
% num_rows = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tot_inds = mbp_inds_ab_array(symm_orders, num_rows);
% ct1 = num_rows;
% coeffs = zeros(num_rows,1)
% while ct1 > 0
%     
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;
% rots1 = rot_mats(:,:,1); r1 = rots1(:,1:3); r2 = rots1(:,4:6);
% M1 = calc_Mfunc_ab_array(symm_orders, num_rows, rots1);
% 
% 
% g_ypi = vrrotvec2mat([0,1,0, pi]);
% r1_ges = g_ypi*r2; r2_ges = g_ypi*r1;
% rots_ges = [r1_ges, r2_ges];
% M1_ges = calc_Mfunc_ab_array(symm_orders, num_rows, rots_ges);
% 
% 
% % nrots = size(rot_mats,3);
% 
% coeffs = col1*(rand(1749,1)*100);
% 
% s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;
% rots1 = rot_mats(:,:,1); r1 = rots1(:,1:3); r2 = rots1(:,4:6);
% M1 = calc_Mfunc_ab_array(symm_orders, num_rows, rots1);
% nrots = size(rot_mats,3);
% 
% g_ypi = vrrotvec2mat([0,1,0, pi]);
% r1_ges = g_ypi*r2; r2_ges = g_ypi*r1;
% rots_ges = [r1_ges, r2_ges];
% M1_ges = calc_Mfunc_ab_array(symm_orders, num_rows, rots_ges);

