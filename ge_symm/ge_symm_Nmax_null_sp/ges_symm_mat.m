clear all; clc;

curr_pwd = split(pwd,'/');
top_dir = '';
for ct1=1:length(curr_pwd)
    top_dir = strcat(top_dir,curr_pwd{ct1},'/');
    if (strcmp(curr_pwd{ct1},'MATLAB_Codes'))
        break;
    end
end
util_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(util_dir));

s1 = set_vars();
Nmax = s1.Nmax; pt_grp = s1.pt_grp;
fname = [top_dir,'data_files', '/ptgrp_',pt_grp];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mat_name = [fname,'/cryst_symm/symm_ab_',...
pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
symm_orders = s1.symm_orders;
nsymm = size(symm_orders,1);
a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1, b1);
num_rows = sum((2*a1+1).*(2*b1+1).*(2*c1+1));

tot_inds = mbp_inds_ab_array(symm_orders, num_rows);


% ges_mat = sparse(num_rows, num_rows);
ges_mat = zeros(num_rows, num_rows);

for ct1=1:num_rows
%     ct1
    a1 = tot_inds(ct1,3);
    b1 = tot_inds(ct1,4);
    gamma1 = tot_inds(ct1,5);
    alpha1 = tot_inds(ct1,6);
    beta1 = tot_inds(ct1,7);
    
    a2 = b1;
    b2 = a1;
    gamma2 = gamma1;
    alpha2 = beta1;
    beta2 = alpha1;
    
    ind1 = find((tot_inds(:,3) == a2) & ...
        (tot_inds(:,4) == b2) & ...
        (tot_inds(:,5) == gamma2) & ...
        (tot_inds(:,6) == alpha2) & ...
        (tot_inds(:,7) == beta2));
    
    ges_mat(ct1,ind1) = (-1)^(a1+b1);
end
ges_mat = ges_mat' - eye(num_rows, num_rows);

col1 = null(ges_mat);

mat_name = [fname,'/ge_symm_null/Y_ges_Nmax_',...
    num2str(Nmax),'.mat'];
save(mat_name, 'col1');


Q_mat = col1*col1';
[v, d] = eig(full(Q_mat));
col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
S = orth(v(:,col));
if any(col)
    S = orth(v(:,col));
end
mat_name = [fname,'/ge_symm_null/Sarr_ges_Nmax_',...
    num2str(Nmax),'.mat'];
save(mat_name,'S');

rmpath(genpath(util_dir));