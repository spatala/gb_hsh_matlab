function S = compute_basis(ga_s,gb_s,n_gen, a_val,b_val, TOL)
gs1 = ga_s{n_gen}; gs2 = gb_s{n_gen};
U_a = rotation_mat(a_val, gs1); U_b = rotation_mat(b_val, gs2);
R1 = sparse((kron(U_a, U_b)).');
S = sp_orth(sp_null(R1 - eye(size(R1)), 1, TOL));
end

