% function [] = check_combining_Nc(top_dir, pt_grp, Nmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 
%%%% 

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
data_fname0 = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ga_s, gb_s, num_gen] = get_symmgen_mats(pt_grp);
symm_orders = zeros(Nmax^2,2);

a_val = 4; b_val = 2; c_val = min(a_val, b_val);

n_gen = 2;


gs1 = ga_s{n_gen}; gs2 = gb_s{n_gen};
Na = 2*a_val; Nb = 2*b_val; Nc = 2*c_val;

R1 = full(transpose(so4_irrep(gs1,gs2,Na,Nb)));
[v1, d1] = eig(R1);
col1 = (abs(imag(diag(d1)))<1e-5 & abs(real(diag(d1))-1)<1e-5);
if any(col1)
    X1 = orth(v1(:,col1));
    X1 = kron(eye(Nc+1), X1);
end

R2 = generate_csymm_ab(gs1, gs2, a_val, b_val);
[v2, d2] = eig(R2);
col2 = (abs(imag(diag(d2)))<1e-5 & abs(real(diag(d2))-1)<1e-5);
if any(col2)
    X2 = orth(v2(:,col2));
end

% nsz = size(X2,2);
% diff_vecs = zeros(nsz,1);
% tA1 = X1;
% for ct1 = 1:nsz
%     tB1 = X2(:,ct1);
%     diff_vecs(ct1) = norm(tA1*(tA1\tB1) - tB1);
% end
