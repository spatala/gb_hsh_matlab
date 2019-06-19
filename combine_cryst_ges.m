clear all; clc;


s1 = load('data_files/ptgrp_C2/cryst_symm/Sarr_4_4.mat');
% X0 = sparse(s1.Sarr_nmax);
% X0 = s1.Sarr_nmax;
X0 = s1.S;
P0 = X0*X0';

s1 = load('ge_symm/Y_ges.mat');
% Y1 = sparse(s1.S);
Y1 = s1.col1;
Y1 = orth(Y1);
Q1 = Y1*Y1';

[v, d] = eig(P0*Q1*P0);
col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
            
if any(col)
    % display(v(:,col));
    S = orth(v(:,col));
    mat_name = 'Sarr_C2_4_4_combined.mat';
    save(mat_name,'S');
end
