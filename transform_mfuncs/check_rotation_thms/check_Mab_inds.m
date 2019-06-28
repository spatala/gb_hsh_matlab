clear all; clc;

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

% s1 = set_vars();
% Nmax = s1.Nmax; pt_grp = s1.pt_grp;
%
% data_fname0 = [top_dir,'data_files/ptgrp_',pt_grp,'/cryst_symm/'];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% num_cols = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% The variable tot_inds contains the mapping between the column-index
% %%%% for the M-function and the (a,b,gamma,alpha, beta) values.
% tot_inds = mbp_inds_ab_array(symm_orders, num_cols);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;
rots1 = rot_mats(:,:,1);
% rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
g1 = rots1(:,1:3); g2 = rots1(:,4:6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mvec_12_nmax = zeros(1,num_cols); Mvec_21_nmax = zeros(1,num_cols);
% ind_start = 1;
% for a_val=0:Nmax
%     for b_val=0:Nmax

a_val = 4; b_val = 3;
c_val = min(a_val, b_val);
Na = 2*a_val; Nb = 2*b_val; Nc = 2*c_val; nsz = (Na+1)*(Nb+1);

R_ab_12 = so4_irrep(g1,g2,Na,Nb);
[Mvec_ab_12, inds, Rab_inds] = calc_Mfunc_Rab_full(R_ab_12, a_val, b_val);
Mvec1_ab_12 = calc_Mfunc_Rab(R_ab_12, a_val, b_val);

gz = vrrotvec2mat([0,0,1,rand()*2*pi]);
Rgz_ab_12 = so4_irrep(gz*g1,gz*g2,Na,Nb);
Mvec2_ab_12 = calc_Mfunc_Rab(Rgz_ab_12, a_val, b_val);

norm(Mvec1_ab_12-Mvec_ab_12)

norm(Mvec1_ab_12-Mvec2_ab_12)


% % %         R_ab_21 = kron(U2, U1);
% % %         Mvec_ab_21 = calc_Mfunc_Rab(R_ab_21, a_val, b_val);
% % % tMvec_ab_12 = calc_Mfunc(a_val, b_val, [g1,g2]);

% % %         Mvec_ab_21 = calc_Mfunc(a_val, b_val, [g2,g1]);
% % norm(Mvec_ab_12 - tMvec_ab_12)
% %
% % % nsz1 = nsz*(Nc+1);
% % % ind_stop = ind_start + nsz1 - 1;
% % % [ind_start, ind_stop]
% % % Mvec_12_nmax(ind_start:ind_stop) = Mvec_ab_12;
% % % Mvec_21_nmax(ind_start:ind_stop) = Mvec_ab_21;
% % % ind_start = ind_stop + 1;
% % %     end
% % % end

rmpath(genpath(util_dir));