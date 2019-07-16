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

mat_name = [data_fname0,'symm_mat_full_ges_nmax_',num2str(Nmax),'.mat'];
% save(mat_name,'symm_mat','symm_orders');

s1 = load(mat_name); symm_mat = s1.symm_mat;

% tic;
% [v,d] = eig(full(symm_mat));
% toc;

% col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
% if any(col)
%     Xarr = orth(v(:,col));
% end

tic;
nsz = size(symm_mat,1);
symm_mat1 = symm_mat - speye(nsz,nsz);
S = spnull(symm_mat1);


% save('Sarr_ges_nmax_4.mat','S');
% s1 = load('')

mat_name = [data_fname0, 'Sarr_abc_combined_csymm_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
X0 = sparse(s1.S);

ne_max = size(X0,2);

P0 = X0*X0';
Q1 = S*S';

R1 = P0*Q1*P0;
R2 = R1 - speye(nsz,nsz);
S = spnull(R2);
toc;

% sig_val = 1.0-1e-6;
% [v,d] = eigs(R1, ne_max, sig_val);

% 
% tic;
% [v,d] = eigs(symm_mat,1000,'largestreal');
% toc;

% tic;
% [Q,R] = qr(symm_mat);
% toc;



display(size(S))
