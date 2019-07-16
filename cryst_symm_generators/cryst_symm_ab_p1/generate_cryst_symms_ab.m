function symm_mat = generate_cryst_symms_ab(top_dir, pt_grp, Nmax, nsymm)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ga_s, gb_s, ~, ~] = get_symmgen_mats(pt_grp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_val = symm_orders(nsymm,1); b_val = symm_orders(nsymm,2);
gs1 = ga_s{nsymm}; gs2 = gb_s{nsymm};
symm_mat = generate_csymm_ab(gs1, gs2, a_val, b_val);

% display(nnz(symm_mat));
% [unique(real(symm_mat(:)));unique(imag(symm_mat(:)))]
% save_symm_mat(symm_orders, symm_mat, data_fname0, ct4);

end

function save_symm_mat(symm_orders, symm_mat, fname, ngen)
Nmax = max(symm_orders(:)); 
mat_name = [fname,'symm_mat_full_ngen_',num2str(ngen)...
    ,'_nmax_',num2str(Nmax),'.mat'];
save(mat_name,'symm_mat','symm_orders');
end





