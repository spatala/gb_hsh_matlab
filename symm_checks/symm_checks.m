function [] = symm_checks(pt_grp, Nmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Code to check the rotation of SO(4) Functions R^(a,b)
%%%%%%%
% clear all; clc;
% pt_grp = 'C1'; Nmax = 1;
% pt_grp = 'Oh'; Nmax = 4;

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
% [symm_rots, ~] = get_symm_rots(g1,g2, pt_grp, data_fname,1);
[symm_rots, ~] = get_symm_rots(g1,g2, pt_grp, data_fname,2);


nrots = size(symm_rots,3);
q1 = zeros(nrots,1); Q1  = zeros(nrots,3);
a1 = zeros(nrots,1); lb1 = zeros(nrots,1);
q2 = zeros(nrots,1); Q2  = zeros(nrots,3);
a2 = zeros(nrots,1); lb2 = zeros(nrots,1);

for ct1 = 1:nrots
    g1_1 = symm_rots(:,1:3,ct1);  axang1 = vrrotmat2vec(g1_1');
    ax1 = axang1(1:3); w_1 = axang1(4);
    
    q1(ct1,1) = cos(w_1/2); Q1(ct1,:) = ax1*sin(w_1/2);
    a1(ct1,1) = q1(ct1,1)-1i*Q1(ct1,3); lb1(ct1,1) = -Q1(ct1,2)-1i*Q1(ct1,1);
    
    g1_2 = symm_rots(:,4:6,ct1); axang2 = vrrotmat2vec(g1_2');
    ax2 = axang2(1:3); w_2 = axang2(4);

    q2(ct1,1) = cos(w_2/2); Q2(ct1,:) = ax2*sin(w_2/2);
    a2(ct1,1) = q2(ct1,1)-1i*Q2(ct1,3); lb2(ct1,1) = -Q2(ct1,2)-1i*Q2(ct1,1);
end

s = struct();
s.q1 = q1; s.Q1 = Q1; s.a1 = a1; s.lb1 = lb1;
s.q2 = q2; s.Q2 = Q2; s.a2 = a2; s.lb2 = lb2;

mat_name = [data_fname0,'Sarr_MabInds_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);

tot_Uprops = s1.tot_Uprops;
vec_inds = tot_inds(unique(tot_Uprops(:,1)),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_gbnull_nmax_',num2str(Nmax),'.mat'];
% mat_name = [data_fname0, ...
%     'Sarr_cryst_ges_nmax_',num2str(Nmax),'.mat'];
% mat_name = [data_fname0,...
%     'Sarr_abc_combined_csymm_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
nsymm_evs = size(S,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


num_rows = size(S,1);
tic;
Mvec = compute_Mvec(s,nrots, num_rows, vec_inds);
toc;

Smvec = Mvec*S;
disp(norm(full(Smvec - Smvec(1,:)))/size(S,2))

% disp(size(Mvec))

% g1_1 = symm_rots(:,1:3,1); g1_2 = symm_rots(:,4:6,1);
% ma1 = rots_to_angs(g1_1, g1_2);
% 
% a_val = symm_orders(:,1)'; b_val = symm_orders(:,2)';
% c_val = min(a_val, b_val);
% num_cols = sum((2*a_val+1).*(2*b_val+1).*(2*c_val+1));
% 
% Mvec1 = zeros(1,num_cols);
% for a=a_val
%     for b=b_val
%         M1 = mbp_basis(a, b, [ma1(1), ma1(2), ma1(3), ma1(4), ma1(5)]);
%         
%         cond1 = tot_inds(:,3)==a & tot_inds(:,4)==b;
%         ind_start = find(cond1,1);
%         ind_stop  = find(cond1,1,'last');
%         
%         Mvec1(ind_start:ind_stop) = M1;
%     end
% end
% 
% nsymm_rots = size(symm_rots,3);
% diff_vec = zeros(nsymm_evs*(nsymm_rots-1),1);
% ct3 = 1;
% for ct2=2:nsymm_rots
%     ct2
%     g2_1 = symm_rots(:,1:3,ct2); g2_2 = symm_rots(:,4:6,ct2);
%     ma2 = rots_to_angs(g2_1, g2_2);
%     for ct1=1:nsymm_evs
%         Mvec2 = zeros(1,num_cols);
%         for a=a_val
%             for b=b_val
%                 M2 = mbp_basis(a, b, [ma2(1), ma2(2), ma2(3), ma2(4), ma2(5)]);       
%                 cond1 = tot_inds(:,3)==a & tot_inds(:,4)==b;
%                 ind_start = find(cond1,1); ind_stop  = find(cond1,1,'last');
%                 Mvec2(ind_start:ind_stop) = M2;
%             end
%         end
%         if norm(Mvec1*S(:,ct1) - Mvec2*S(:,ct1)) > 1e-12
%             disp([ct2, ct1])
%         end
%         diff_vec(ct3) = norm(Mvec1*S(:,ct1) - Mvec2*S(:,ct1)); ct3 = ct3 + 1;
%     end
% end
% disp(max(abs(diff_vec)));

rmpath(genpath(util_dir));
end