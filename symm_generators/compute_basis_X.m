function S = compute_basis_X(S0, ga_s,gb_s,n_gen, a_val,b_val, TOL)
gs1 = ga_s{n_gen}; gs2 = gb_s{n_gen};
U_a = rotation_mat(a_val, gs1); U_b = rotation_mat(b_val, gs2);
R1 = sparse((kron(U_a, U_b)).');
S = S0 * sp_orth(sp_null(S0' * R1 * S0 - eye(size(S0, 2)), 1, TOL));
end
