function [] = generate_cryst_symm_ab(top_dir, pt_grp, Nmax, TOL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'aPLUSb_max_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ga_s, gb_s, num_gen, Laue] = get_symmgen_angs(pt_grp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_inds = ((Nmax+1)*(Nmax+2)/2);

symm_orders = zeros(n_inds,2);
ab_inds = zeros(n_inds,2);
tct3 = 1;
for tct1=0:Nmax
    for tct2=0:tct1
        ab_inds(tct3,:) = [tct2,tct1-tct2];
        tct3 = tct3 + 1;
    end
end

Sab_arr = {};
num_rows = 0; num_cols = 0; ns_ord = 1;

for tct1 = 1:n_inds
    a_val = ab_inds(tct1,1); b_val = ab_inds(tct1,2);
    disp([a_val, b_val])
    for ct4 = 1:2*num_gen
        if ct4 == 1
            S0 = compute_basis(ga_s,gb_s,ct4,a_val,b_val, TOL);
        else
            S0 = compute_basis_X(S0, ga_s,gb_s,ct4,a_val,b_val, TOL);
        end
        if size(S0,2) == 0, break; end
    end
    c_val = min(a_val, b_val); Nc = 2*c_val;
    S0 = kron(eye(Nc+1), S0);
    if Laue
        R1 = generate_ypi_left_ab(a_val,b_val);
        S0 = S0 * sp_orth(sp_null(S0' * R1 * S0 - eye(size(S0, 2)), 1, TOL));
    end
    
    if (size(S0,2) > 0)
        Na = 2*a_val; Nb = 2*b_val;
        num_rows = num_rows + (Na+1)*(Nb+1)*(Nc+1);
        num_cols = num_cols + size(S0,2);
        Sab_arr{ns_ord} = S0;
        symm_orders(ns_ord,:) = [a_val, b_val]; ns_ord = ns_ord + 1;
    end
end
symm_orders(ns_ord:end,:) = [];
nsymm = ns_ord-1;

S1_arr = sparse(zeros(num_rows, num_cols));
row_ind_start = 1;
col_ind_start = 1;
for ns_ord = 1:nsymm    
    Sab = Sab_arr{ns_ord};
    row_ind_stop = row_ind_start + size(Sab,1) - 1;
    col_ind_stop = col_ind_start + size(Sab,2) - 1;    
    S1_arr(row_ind_start:row_ind_stop, col_ind_start:col_ind_stop) = Sab;
    row_ind_start = row_ind_stop+1;
    col_ind_start = col_ind_stop+1;
end
S = clean(S1_arr, TOL);
mat_name = [data_fname0,'Sarr_abc_combined_csymm_aPLUSb_max_',num2str(Nmax),'.mat'];
save(mat_name, 'S','-v7.3');

mat_name = [data_fname0,'symm_ab_',pt_grp,'_aPLUSb_max_',num2str(Nmax),'.mat'];
save(mat_name,'symm_orders','TOL','-v7.3');

end

function S = compute_basis(ga_s, gb_s, n_gen, a_val, b_val, TOL)
gs1 = ga_s{n_gen}; gs2 = gb_s{n_gen};
U_a = rotation_mat(a_val, gs1); U_b = rotation_mat(b_val, gs2);
R1 = sparse((kron(U_a, U_b)).');
S = sp_orth(sp_null(R1 - eye(size(R1)), 1, TOL));
end


function S = compute_basis_X(S0, ga_s,gb_s,n_gen, a_val,b_val, TOL)
gs1 = ga_s{n_gen}; gs2 = gb_s{n_gen};
U_a = rotation_mat(a_val, gs1); U_b = rotation_mat(b_val, gs2);
R1 = sparse((kron(U_a, U_b)).');
S = S0 * sp_orth(sp_null(S0' * R1 * S0 - eye(size(S0, 2)), 1, TOL));
end
