function [] = check_basis_orth(pt_grp, Nmax, coeffs_typ)
% clear all; clc;
% pt_grp = 'Oh'; Nmax = 8; coeffs_typ = 'aPLUSb_max'; coeffs_typ = 'nmax';

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_',coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mat_name = [data_fname0, ...
    'Sarr_cryst_ges_gbnull_nmax_',num2str(Nmax),'.mat'];

s1 = load(mat_name); S = (s1.S);
nsymm_evs = size(S,2);

Sval = zeros(nsymm_evs,nsymm_evs);
for ct1 = 1:nsymm_evs
    Svec1 = full(S(:,ct1));
    for ct2 = 1:nsymm_evs
        Svec2 = full(S(:,ct2));
        Sval(ct1,ct2) = sum(Svec1.*Svec2);
    end
end

norm(Sval - eye(nsymm_evs))/(nsymm_evs)

rmpath(genpath(util_dir));

end

