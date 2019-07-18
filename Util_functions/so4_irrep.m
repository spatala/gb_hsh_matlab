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
%%%%    SO(4) irrep - square matrix with (Na+1)*(Nb+1) rows.
%%%% 

%%%%% If Identity matrix, then no need to call rotation() function!
tm1 = g1 - eye(3); tm1 = tm1(:);
tm2 = g2 - eye(3); tm2 = tm2(:);

if (norm(tm1) < 1e-14)
    U1 = eye(Na+1);
else
    ax_ang_1 = vrrotmat2vec(g1'); 
    U1 = rotation( ax_ang_1(1:3),  ax_ang_1(4), Na);
end

if (norm(tm2) < 1e-14)
    U2 = eye(Nb+1);
else
    ax_ang_2 = vrrotmat2vec(g2');
    U2 = rotation( ax_ang_2(1:3),  ax_ang_2(4), Nb);
end
R_ab_12 = sparse(kron(U1, U2));
end