function [] = combine_cryst_ge_symms(top_dir, pt_grp, Nmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mat_name = [data_fname0,'symm_mat_full_ges_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_mat = s1.symm_mat;

nsz = size(symm_mat,1);
symm_mat1 = symm_mat - speye(nsz,nsz);
Y1 = spnull(symm_mat1);

mat_name = [data_fname0, 'Sarr_abc_combined_csymm_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
X0 = sparse(s1.S);

ne_max = size(X0,2);

P0 = X0*X0';
Q1 = Y1*Y1';

R1 = P0*Q1*P0;
R2 = R1 - speye(nsz,nsz);
S = spnull(R2);
toc;

display(size(Y1))
end