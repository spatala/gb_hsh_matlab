clear all; clc;

curr_pwd = split(pwd,'/');
top_dir = '';
for ct1=1:length(curr_pwd)
    top_dir = strcat(top_dir,curr_pwd{ct1},'/');
    if (strcmp(curr_pwd{ct1},'MATLAB_Codes'))
        break;
    end
end
top_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(top_dir));

s1 = set_vars();
Nmax = s1.Nmax;
pt_grp = s1.pt_grp;


rmpath(genpath(top_dir));
% s1 = load('../data_files/ptgrp_C2/cryst_symm/Sarr_8_8.mat');
% X0 = s1.S;
% P0 = X0*X0';
% 
% 
% s1 = load('../ge_symm/Y_ges.mat');
% Y1 = s1.col1;
% % Y1 = orth(Y1);
% Q1 = Y1*Y1';
% % 
% [v, d] = eig(P0*Q1*P0);
% col = (abs(imag(diag(d)))<1e-5 & abs(real(diag(d))-1)<1e-5);
%             
% if any(col)
%     % display(v(:,col));
%     S = orth(v(:,col));
%     mat_name = 'Sarr_C2_4_4_combined.mat';
%     save(mat_name,'S');
% end
