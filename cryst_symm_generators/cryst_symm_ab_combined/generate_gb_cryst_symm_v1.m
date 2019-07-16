function [] = generate_gb_cryst_symm_v1()
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

s1 = set_vars(); Nmax = s1.Nmax; pt_grp = s1.pt_grp;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, ~, num_gen, Laue] = get_symmgen_mats(pt_grp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
nsymm = size(symm_orders,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Laue
    num_symm = 2*num_gen + 1;
else
    num_symm = 2*num_gen;
end

Sab_arr = cell(nsymm,1);

num_rows = 0; num_cols = 0;
% nproj = num_symm;
for ns_ord = 1:nsymm
    X1 = cell(num_symm,1);
    a_val = symm_orders(ns_ord,1); b_val = symm_orders(ns_ord,2);
    c_val = min(a_val, b_val); Nc = 2*c_val;
    
    ct1 = 1;
    mat_name = ['Sarr_ab_',num2str(a_val),'_',num2str(b_val),'_',num2str(ct1),'.mat'];
    s1 = load(mat_name); X1{ct1} = s1.S;
    if Laue
        for ct1 = 2:num_symm-1
            ct0 = ct1-1;
            mat_name = ['Sarr_ab_',num2str(a_val),'_',num2str(b_val),'_',num2str(ct1),'.mat'];
            X1{ct1} = get_Xarr_proj(X1{ct0}, mat_name);
        end
        X1{ct1} = kron(eye(Nc+1), X1{ct1});
        ct1 = ct1 + 1; ct0 = ct1-1;
        mat_name = ['Sarr_ab_',num2str(a_val),'_',num2str(b_val),'_',num2str(ct1),'.mat'];
        X1{ct1} = get_Xarr_proj(X1{ct0}, mat_name); 
    else
        for ct1 = 2:num_symm
            ct0 = ct1-1;
            mat_name = ['Sarr_ab_',num2str(a_val),'_',num2str(b_val),'_',num2str(ct1),'.mat'];
            X1{ct1} = get_Xarr_proj(X1{ct0}, mat_name);
        end
        X1{ct1} = kron(eye(Nc+1), X1{ct1});
    end
    
    Sab = X1{ct1};
    
    Na = 2*a_val; Nb = 2*b_val;
    c_val = min(a_val, b_val); Nc = 2*c_val;
    num_rows = num_rows + (Na+1)*(Nb+1)*(Nc+1);
    num_cols = num_cols + size(Sab,2);
    
    Sab_arr{ns_ord} = Sab;
end



S1_arr = zeros(num_rows, num_cols);
row_ind_start = 1;
col_ind_start = 1;
for ns_ord = 1:nsymm
    a_val = symm_orders(ns_ord,1); b_val = symm_orders(ns_ord,2);
%     Na = 2*a_val; Nb = 2*b_val;
%     c_val = min(a_val, b_val); Nc = 2*c_val;
    
    Sab = Sab_arr{ns_ord};
    
%     Sabc = kron(eye(Nc+1), Sab);
    
    
    row_ind_stop = row_ind_start + size(Sab,1) - 1;
    col_ind_stop = col_ind_start + size(Sab,2) - 1;
    
    S1_arr(row_ind_start:row_ind_stop, col_ind_start:col_ind_stop) = Sab;
    
    row_ind_start = row_ind_stop+1;
    col_ind_start = col_ind_stop+1;
end

S1_arr = clean(S1_arr);

nproj = 5;
s1 = load(['Sarr_cryst_combined_abc_nmax_4_',num2str(nproj),'.mat']);
S2 = s1.S;

% % S1_arr = orth(S1_arr);
% check_equi_basis(S1_arr,S2)
% check_equi_basis(S2,S1_arr)
%
% X0 = clean(S1_arr);
%
% fname1 = [top_dir,'data_files/ptgrp_Oh/nmax_4'];
% fname2 = '/symm_mat_full_ngen_5_nmax_4.mat';
% mat_name = [fname1, fname2];
% s1 = load(mat_name);
% symm_mat = s1.symm_mat;
% [v,d] = eig(full(symm_mat));
% col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
% if any(col)
%     Y1 = orth(v(:,col));
%     Y1 = clean(Y1);
% end
% P0 = X0*X0';
% Q1 = Y1*Y1';
% R1 = P0*Q1*P0;
% [v, d] = eig(R1);
%
% col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
% if any(col)
%     X0 = orth(v(:,col));
%     X0 = clean(X0);
% end
%
% fname2 = '/symm_mat_full_ngen_6_nmax_4.mat';
% mat_name = [fname1, fname2];
% s1 = load(mat_name);
% symm_mat = s1.symm_mat;
% [v,d] = eig(full(symm_mat));
% col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
% if any(col)
%     Y1 = orth(v(:,col));
%     Y1 = clean(Y1);
% end
% P0 = X0*X0';
% Q1 = Y1*Y1';
% R1 = P0*Q1*P0;
% [v, d] = eig(R1);
%
% col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
% if any(col)
%     X0 = orth(v(:,col));
%     X0 = clean(X0);
% end

rmpath(genpath(util_dir));

end

function Xarr = get_Xarr_proj(X1, mat_name)
P1 = X1*X1';
s1 = load(mat_name); Y1 = s1.S;
Q1 = Y1*Y1';
R1 = P1*Q1*P1;
[v,d] = eig(R1);
col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
if any(col)
    Xarr = orth(v(:,col));
end
end
