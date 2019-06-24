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

symm_orders = [3,2];
a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1,b1);
num_rows = sum((2*a1+1).*(2*b1+1).*(2*c1+1));

s1 = load([top_dir,'GB_Parameters/rand_gb_rots.mat']); rot_mats = s1.rot_mats;
rots1 = rot_mats(:,:,1); r1 = rots1(:,1:3); r2 = rots1(:,4:6);
M1 = calc_Mfunc_ab_array(symm_orders, num_rows, rots1);
nrots = size(rot_mats,3);

diff_ct2 = zeros(nrots-1,1);
for ct2 = 2:nrots
    rots2 = rot_mats(:,:,ct2); t1 = rots2(:,1:3); t2 = rots2(:,4:6);
    % r1t1 = r1*(t1'); r2t2 = r2*(t2');
    r1t1 = t1'*r1; r2t2 = t2'*r2;
    M1t = calc_Mfunc_ab_array(symm_orders, num_rows, [r1t1, r2t2]);
    
%     t1_axang = vrrotmat2vec(t1'); t2_axang = vrrotmat2vec(t2');
    t1_axang = vrrotmat2vec(t1'); t2_axang = vrrotmat2vec(t2');
    
    Na = 2*a1; Nb = 2*b1; Nc = 2*c1;
    R1 = kron(rotation(t1_axang(1:3),t1_axang(4),Na),...
        rotation(t2_axang(1:3),t2_axang(4),Nb));
    R1 = kron(eye(Nc+1), R1);
    
    diff_ct2(ct2) = norm(M1*R1 - M1t);
end

% diff_vec = zeros(num_symm*num_symm,1);
% ct3 = 1;
% for ct1 = 1:num_symm
%     for ct2 = 1:num_symm
%         S1 = SymmMat{ct1}; S2 = SymmMat{ct2};
%
%         r1s = r1*S1; r2s = r2*S2;
%         rots_symm = [r1s, r2s];
%         M1_symm = calc_Mfunc_ab_array(symm_orders, num_rows, rots_symm);
%
%         diff_vec(ct3) = norm(M1_symm*Svec - M1*Svec);
%         ct3 = ct3 + 1;
%     end
% end
%
% max(diff_vec)

rmpath(genpath(top_dir));