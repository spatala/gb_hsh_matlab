function [] = generate_ge_symms(top_dir, pt_grp, Nmax)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ypi_left = generate_ypi_left(symm_orders);
flip_mat = generate_flip_mat(symm_orders);
symm_mat = ypi_left*flip_mat;

% Nmax = max(symm_orders(:)); 
% mat_name = [data_fname0,'symm_mat_full_ges_nmax_',num2str(Nmax),'.mat'];
% save(mat_name,'symm_mat','symm_orders');

mat_name = [data_fname0, ...
    'Sarr_ges_nmax_',num2str(Nmax),'.mat'];
nsz = size(symm_mat,1);
symm_mat1 = symm_mat - speye(nsz,nsz);
S = spnull(symm_mat1);
save(mat_name,'S');


end




