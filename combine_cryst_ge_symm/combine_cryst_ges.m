clear all; clc;

addpath(genpath('../Util_functions/'));
addpath(genpath('../GB_Parameters/'));

fname = get_dir_name();
% pt_grp = 'C2';
pt_grp = 'O';
Nmax = 8;

mat_name = [fname,'/ptgrp_',pt_grp,'/cryst_symm/symm_ab_',...
    pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
symm_orders = s1.symm_orders;
nsymm = size(symm_orders,1);
a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1, b1);
num_rows = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
Nmax = max(a1);

s1 = load([fname,'/ptgrp_',pt_grp,'/cryst_symm/Sarr_Nmax_',...
    num2str(Nmax),'.mat']);
X0 = s1.S;
for ct1=1:size(X0,2)
    ind1 = find(sqrt(real(X0(:,ct1)).^2 + imag(X0(:,ct1)).^2) < 1e-10);
    X0(ind1,ct1) = 0;
end
X0 = sparse(X0);
P0 = X0*X0';

% s1 = load([fname,'/ptgrp_',pt_grp,'/ge_symm/Y_ges_Nmax_',...
%     num2str(Nmax),'.mat']);
% Y1 = sparse(s1.col1);
% Y1 = orth(Y1);

s1 = load([fname,'/ptgrp_',pt_grp,'/ge_symm/Sarr_ges_Nmax_',...
    num2str(Nmax),'.mat']);
Y1 = s1.S;
Q1 = Y1*Y1';

eig_mat = P0*Q1*P0;

[v, d] = eig(full(eig_mat));
col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
            
if any(col)
    S = orth(v(:,col));
    mat_name = [fname,'/ptgrp_',pt_grp,'/combined_symm/Sarr_combined_Nmax_',...
    num2str(Nmax),'.mat'];
    save(mat_name,'S');
end

rmpath(genpath('../Util_functions/'));
rmpath(genpath('../GB_Parameters/'));
