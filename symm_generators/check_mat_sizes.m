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

%%%%
gen_symm_orders(top_dir, pt_grp, Nmax);
%%%%

data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%

% for ngen = 1:6
for ngen = [2,4]
    ngen
    mat_name = [data_fname0,'symm_mat_full_ngen_',num2str(ngen)...
        ,'_nmax_',num2str(Nmax),'.mat'];
%     mat_name = [data_fname0,'symm_mat_full_ngen_2_nmax_',num2str(Nmax),'.mat'];
    display(mat_name) 
    s1 = load(mat_name);
    R1 = full(s1.symm_mat);
    
    % num_eps = 8;
    % num_nz_elem = zeros(num_eps,1);
    % for ct1=1:num_eps
    %     eps_tol = (1/(10^(ct1+8)));
    %     num_nz_elem(ct1) = size(R1(:),1) - size(find(abs(R1(:)) < eps_tol),1);
    % end
    % display(num_nz_elem)
    
    [v0, d0] = eig(R1);
    
    col0 = (abs(imag(diag(d0)))<1e-5 & abs(real(diag(d0))-1)<1e-5);
    if any(col0)
        Y0 = orth(v0(:,col0));
%         Q0 = Y0*Y0';
    end
    
%     Y0 = clean(Y0);
    
    num_eps = 8;
    num_nz_elem = zeros(num_eps,1);
    for ct1=1:num_eps
        eps_tol = (1/(10^(ct1+8)));
    
        num_nz_elem(ct1) = size(Y0(:),1) - size(find(abs(Y0(:)) < eps_tol),1);
    end
    display(size(Y0,2))
    
%     num_eps = 8;
%     num_nz_elem = zeros(num_eps,1);
%     for ct1=1:num_eps
%         eps_tol = (1/(10^(ct1+8)));
%         num_nz_elem(ct1) = size(Q0(:),1) - size(find(abs(Q0(:)) < eps_tol),1);
%     end
    
    display(num_nz_elem)
    display('+++++++++++++++++++++')
    % num_nz_elem;
    
end
%
% [v,d] = eig(R1);
