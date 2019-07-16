function [] = generate_gb_cryst_symm(pt_grp, Nmax)
% clear all; clc;

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

% s1 = set_vars(); Nmax = s1.Nmax; pt_grp = s1.pt_grp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%
% gen_symm_orders(top_dir, pt_grp, Nmax);
%%%%
% generate_ab_symms(top_dir, pt_grp, Nmax);
%%%%
% combine_cryst_symm_ab(top_dir, pt_grp, Nmax);
%%%%
% generate_ge_symms(top_dir, pt_grp, Nmax);
%%%%
% combine_cryst_ges(data_fname0, Nmax)
%%%%

rmpath(genpath(util_dir));

end

function combine_cryst_ges(data_fname0, Nmax)
mat_name = [data_fname0, ...
    'Sarr_ges_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); Y1 = s1.S;

nsz = size(Y1,1);

mat_name = [data_fname0, ...
    'Sarr_abc_combined_csymm_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); X0 = sparse(s1.S);

P0 = X0*X0';
Q1 = Y1*Y1';

R1 = P0*Q1*P0;
R2 = R1 - speye(nsz,nsz);
S = spnull(R2);

mat_name = [data_fname0, ...
    'Sarr_cryst_ges_nmax_',num2str(Nmax),'.mat'];
save(mat_name,'S');

end
