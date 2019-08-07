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

pt_grp = 'Oh'; Nmax = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% m1 = vrrotvec2mat([1,1,1, pi/3]);
m1 = vrrotvec2mat([1,0,0, pi/4]);
% m1 = vrrotvec2mat([1,0,0, 0]);
% q1 = [cos(pi/4), 0, 0, sin(pi/4)];
% m1 = quat2mat(q1);
num = 250;
[X, Y, Z, rots] = gen_gb_bpl_rots(num, m1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
mat_name = [data_fname0, 'Sarr_cryst_ges_gbnull_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
tot_inds = mbp_inds_ab_array(symm_orders);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nrots = size(rots,3);

a_val = symm_orders(:,1)'; b_val = symm_orders(:,2)';
c_val = min(a_val, b_val);
num_cols = sum((2*a_val+1).*(2*b_val+1).*(2*c_val+1));

q1 = zeros(nrots,1); Q1  = zeros(nrots,3);
a1 = zeros(nrots,1); lb1 = zeros(nrots,1);
q2 = zeros(nrots,1); Q2  = zeros(nrots,3);
a2 = zeros(nrots,1); lb2 = zeros(nrots,1);

for ct1 = 1:nrots
    g1_1 = rots(:,1:3,ct1);  axang1 = vrrotmat2vec(g1_1');
    ax1 = axang1(1:3); w_1 = axang1(4);
    
    q1(ct1,1) = cos(w_1/2); Q1(ct1,:) = ax1*sin(w_1/2);
    a1(ct1,1) = q1(ct1,1)-1i*Q1(ct1,3); lb1(ct1,1) = -Q1(ct1,2)-1i*Q1(ct1,1);
    
    g1_2 = rots(:,4:6,ct1); axang2 = vrrotmat2vec(g1_2');
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

num_rows = size(S,1);
% tic;
Mvec = compute_Mvec_v1(s,nrots, num_rows, vec_inds);
% toc;

mat_name = ['Mvec_',pt_grp,'_Nmax_',num2str(Nmax),'_num_',num2str(num),'.mat'];
save(mat_name, 'Mvec','m1','-v7.3');

% Mdiff = (Mvecs(:,vec_inds(:,1)) - Mvec2(:,vec_inds(:,1)));
% 
% num_v = 0;
% for ct1 = 1:nrots
%     if norm(Mdiff(ct1,:))>1e-10
%         num_v = num_v + 1;
%     end
% end
% a1_rots  = s.a1;  a2_rots  = s.a2;
% lb1_rots = s.lb1; lb2_rots = s.lb2;
% q1_rots  = s.q1;  q2_rots  = s.q2;
% Q1_rots  = s.Q1;  Q2_rots  = s.Q2;


% Mvecs1 = sparse(Mvecs);

% save(mat_name, 'Mvecs1', '-v7.3');
% symm_Mvec = Mvecs*S(:,symm_num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure('Position',[1200,50,1200,1200]); hold on;
% fval = abs(symm_Mvec);
% fval = reshape(fval, [num+1, num+1]);
% surf(X,Y,Z, fval);
% shading interp;
% axis equal; axis off;
% % view([1,0,0])
% view([1,1,1]);


rmpath(genpath(util_dir));