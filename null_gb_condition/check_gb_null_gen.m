function [] = check_gb_null_gen(pt_grp, Nmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tot_inds = mbp_inds_ab_array(symm_orders);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_name = [data_fname0, ...
%     'Sarr_gbnull_nmax_',num2str(Nmax),'.mat'];
% mat_name = [data_fname0, ...
%     'Sarr_cryst_ges_gbnull_nmax_',num2str(Nmax),'.mat'];
% 
% s1 = load(mat_name); S = (s1.S);


null_mat = generate_gbnull_mat(top_dir,pt_grp,Nmax);
S = sp_null(null_mat,1,1e-14);

nsymm_evs = size(S,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_vec = zeros(nsymm_evs,1);

w_m  = 0; th_m = pi * rand(); ph_m = 2. * pi * rand();
axang_m = [sin(th_m)*cos(ph_m), sin(th_m)*sin(ph_m), cos(th_m), w_m];
gm = vrrotvec2mat(axang_m);

w_b  = 8 * pi * (rand() - 1); th_b = pi / 2; ph_b = 2. * pi * rand();
axang_b = [sin(th_b)*cos(ph_b), sin(th_b)*sin(ph_b), cos(th_b), w_b];
gb = vrrotvec2mat(axang_b);

a_val = symm_orders(:,1)'; b_val = symm_orders(:,2)';
c_val = min(a_val, b_val);
num_cols = sum((2*a_val+1).*(2*b_val+1).*(2*c_val+1));

Mvec = zeros(1,num_cols);

for ct_ab = 1:size(symm_orders,1)
    a = a_val(ct_ab); b = b_val(ct_ab);
    M1 = mbp_basis(a, b, [w_m, th_m, ph_m, w_b, ph_b]);
    
    cond1 = tot_inds(:,3)==a & tot_inds(:,4)==b;
    ind_start = find(cond1,1);
    ind_stop  = find(cond1,1,'last');
    
    Mvec(ind_start:ind_stop) = M1;
end


for ct1=1:nsymm_evs 
    diff_vec(ct1) = norm(Mvec*S(:,ct1));
end
max(abs(diff_vec))

rmpath(genpath(util_dir));