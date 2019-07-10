function R_ab_12 = so4_irrep(g1,g2,Na,Nb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Function to compute SO(4) Irreducible representative with g1, g2
%%%% rotation matrices
%%%% 
%%%% Input:
%%%% g1, g2:
%%%%    3X3 rotation matrices.
%%%% 
%%%% Na, Nb:
%%%%    Na+1, Nb+1 are the sizes of the SO(3) irreducible representatives.
%%%% 
%%%% Output:
%%%% R_ab_12:
%%%%    SO(4) irrep - (Na+1)*(Nb+1) square matrix 
%%%% 
ax_ang_1 = vrrotmat2vec(g1'); ax_ang_2 = vrrotmat2vec(g2');
U1 = rotation( ax_ang_1(1:3),  ax_ang_1(4), Na);
% U1 = adjust_complex_mats(U1);
U2 = rotation( ax_ang_2(1:3),  ax_ang_2(4), Nb);
% U2 = adjust_complex_mats(U2);
R_ab_12 = kron(U1, U2);
% R_ab_12 = adjust_complex_mats(R_ab_12);
end