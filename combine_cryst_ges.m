clear all; clc;

s1 = load('cryst_symm/Sarr_Nmax_4_C2.mat');
% X0 = sparse(s1.Sarr_nmax);
X0 = s1.Sarr_nmax;
P0 = X0*X0';

s1 = load('ge_symm/Y_ges.mat');
% Y1 = sparse(s1.S);
Y1 = s1.col1;
Y1 = orth(Y1);
Q1 = Y1*Y1';

[v, d] = eig(P0*Q1*P0);
col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
            
if any(col)
    display(v(:,col))
end
