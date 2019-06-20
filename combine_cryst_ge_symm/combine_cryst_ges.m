clear all; clc;

addpath(genpath('../Util_functions/'));
addpath(genpath('../GB_Parameters/'));

Nmax = 4;
fname = get_dir_name();
pt_grp = 'C2';

s1 = load([fname,'/ptgrp_',pt_grp,'/cryst_symm/Sarr_Nmax_',...
    num2str(Nmax),'.mat']);
X0 = s1.Sarr_nmax;
P0 = X0*X0';

s1 = load([fname,'/ptgrp_',pt_grp,'/ge_symm/Y_ges_Nmax_',...
    num2str(Nmax),'.mat']);
Y1 = s1.col1;
Y1 = orth(Y1);
Q1 = Y1*Y1';

[v, d] = eig(P0*Q1*P0);
col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
            
if any(col)
    % display(v(:,col));
    S = orth(v(:,col));
    mat_name = [fname,'/ptgrp_',pt_grp,'/combined_symm/Sarr_combined_Nmax_',...
    num2str(Nmax),'.mat'];
%     mat_name = 'Sarr_C2_4_4_combined.mat';
    save(mat_name,'S');
end

rmpath(genpath('../Util_functions/'));
rmpath(genpath('../GB_Parameters/'));
