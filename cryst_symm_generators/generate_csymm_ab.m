function csymm_mat = generate_csymm_ab(gs1, gs2, a1, b1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Function to compute the full Symmetry matrix that operates on
%%%% M(0,0)....M(Nmax, Nmax) given by symm_orders
%%%%
%%%% Input
%%%% gs1, gs2:
%%%%    3X3 rotation matrices operating on g1, g2.
%%%%
%%%% symm_orders:
%%%%
%%%% Laue:
%%%%
%%%%
%%%% Output
%%%% csymm_mat:
%%%%    Symmetry Rotation operation on M(0,0)....M(Nmax, Nmax)
%%%%

c1 = min(a1,b1); Na = 2*a1; Nb = 2*b1; Nc = 2*c1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rr_ab_12 = so4_irrep(gs1,gs2,Na,Nb); trRr_ab_12 = transpose(Rr_ab_12);
csymm_mat = kron(eye(Nc+1), trRr_ab_12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

