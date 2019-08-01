%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Code to check the rotation of SO(4) Functions R^(a,b)
%%%%%%%
clear all; clc;

% pt_grp = 'C1'; Nmax = 1;
pt_grp = 'Oh'; Nmax = 6;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = load([top_dir,'data_files/GB_Parameters/rand_gb_rots.mat']); 
rot_mats = s1.rot_mats;
rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
g1 = rots1(:,1:3); g2 = rots1(:,4:6);
% M1 = vrrotvec2mat([1,1,1,pi/3]);
% th = rand()*pi; ph = rand()*2*pi;
% nvec = [sin(th)*cos(ph); sin(th)*sin(ph); cos(th)];
% rots1 = mbp_to_rots([M1, nvec]);
% g1 = rots1(:,1:3); g2 = rots1(:,4:6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
tot_inds = mbp_inds_ab_array(symm_orders);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[symm_rots, Laue] = get_symm_rots(g1,g2, pt_grp, data_fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_gbnull_nmax_',num2str(Nmax),'.mat'];
% mat_name = [data_fname0,...
%     'Sarr_abc_combined_csymm_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
nsymm_evs = size(S,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


g1_1 = symm_rots(:,1:3,1); g1_2 = symm_rots(:,4:6,1);
ma1 = rots_to_angs(g1_1, g1_2);

a_val = symm_orders(:,1)'; b_val = symm_orders(:,2)';
c_val = min(a_val, b_val);
num_cols = sum((2*a_val+1).*(2*b_val+1).*(2*c_val+1));

Mvec1 = zeros(1,num_cols); Mvec2 = zeros(1,num_cols);
for a=a_val
    for b=b_val
        M1 = mbp_basis(a, b, [ma1(1), ma1(2), ma1(3), ma1(4), ma1(5)]);
        
        cond1 = tot_inds(:,3)==a & tot_inds(:,4)==b;
        ind_start = find(cond1,1);
        ind_stop  = find(cond1,1,'last');
        
        Mvec1(ind_start:ind_stop) = M1;
    end
end

nsymm_rots = size(symm_rots,3);
diff_vec = zeros(nsymm_evs*(nsymm_rots-1),1);
ct3 = 1;
for ct2=2:nsymm_rots
    ct2
    g2_1 = symm_rots(:,1:3,ct2); g2_2 = symm_rots(:,4:6,ct2);
    ma2 = rots_to_angs(g2_1, g2_2);
    for ct1=1:nsymm_evs
        Mvec2 = zeros(1,num_cols);
        for a=a_val
            for b=b_val
                M2 = mbp_basis(a, b, [ma2(1), ma2(2), ma2(3), ma2(4), ma2(5)]);       
                cond1 = tot_inds(:,3)==a & tot_inds(:,4)==b;
                ind_start = find(cond1,1); ind_stop  = find(cond1,1,'last');
                Mvec2(ind_start:ind_stop) = M2;
            end
        end
        if norm(Mvec1*S(:,ct1) - Mvec2*S(:,ct1)) > 1e-12
            display([ct2, ct1])
        end
        diff_vec(ct3) = norm(Mvec1*S(:,ct1) - Mvec2*S(:,ct1)); ct3 = ct3 + 1;
    end
end

rmpath(genpath(util_dir));