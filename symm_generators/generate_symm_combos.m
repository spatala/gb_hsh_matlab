function [] = generate_symm_combos(top_dir, pt_grp, Nmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 
%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'ge_symm/'];
data_fname1 = [data_fname,'cryst_symm/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Load Grain Exchange Symmetry eigen vectors
mat_name = [data_fname0,'ges_full_symm_orders_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S_ges = s1.S;

%%%% Load Crystal Symmetry eigen vectors
mat_name = [data_fname1,'Sarr_full_symm_orders_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S_csymm = s1.S;

%%%% Combine orthogonal-eigen-vectors
[v1, d1] = combine_XY_symms(S_csymm, S_ges);
col1 = (abs(imag(diag(d1)))<1e-5 & abs(real(diag(d1))-1)<1e-5);

%%%% Save combined (crystal-symmetry and GES) S-vectors
if any(col1)
    save_symm_arr(Nmax, v1, col1, data_fname);
end

end

function [v, d] = combine_XY_symms(X0, Y1)
P0 = X0*X0';
Q1 = Y1*Y1';
[v, d] = eig(P0*Q1*P0);
end

function save_symm_arr(Nmax, v, col, fname)
S = orth(v(:,col));
mat_name = [fname,'Sarr_combined_nmax_',num2str(Nmax),'.mat'];
save(mat_name,'S');
end