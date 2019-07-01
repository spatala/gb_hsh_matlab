function [] = generate_all_symms(top_dir, pt_grp, Nmax)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ga_s, gb_s, num_gen, Laue] = get_symmgen_mats(pt_grp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_symm = 2*num_gen;
for ct4 = 1:num_symm
    ct4
    gs1 = ga_s{ct4}; gs2 = gb_s{ct4};
    symm_mat = generate_csymm(gs1, gs2, symm_orders);
    save_symm_mat(symm_orders, symm_mat, data_fname0, ct4);
end
if (Laue)
    ct4 = ct4 + 1; ct4
    ypi_left = generate_ypi_left(symm_orders);
    save_symm_mat(symm_orders, ypi_left, data_fname0, ct4);    
    
    ct4 = ct4 + 1; ct4
    flip_mat = generate_flip_mat(symm_orders);
    symm_mat = ypi_left*flip_mat;
    save_symm_mat(symm_orders, symm_mat, data_fname0, ct4);
else
    ct4 = ct4 + 1; ct4
    ypi_left = generate_ypi_left(symm_orders);
    flip_mat = generate_flip_mat(symm_orders);
    symm_mat = ypi_left*flip_mat;
    save_symm_mat(symm_orders, symm_mat, data_fname0, ct4);
end
end

function save_symm_mat(symm_orders, symm_mat, fname, ngen)
Nmax = max(symm_orders(:)); 
mat_name = [fname,'symm_mat_full_ngen_',num2str(ngen)...
    ,'_nmax_',num2str(Nmax),'.mat'];
save(mat_name,'symm_mat','symm_orders');
end





