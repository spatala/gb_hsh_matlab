function R_ab_12 = so4_irrep(rot_angs1, rot_angs2, a, b)
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
if (mod(rot_angs1(1),2*pi) < 1e-14)
    U1 = eye(Na+1);
else
    U1 = rotation_mat(a, rot_angs1);    
end

if (mod(rot_angs2(1),2*pi) < 1e-14)
    U2 = eye(Nb+1);
else
    U2 = rotation_mat(b, rot_angs2);
end
R_ab_12 = sparse(kron(U1, U2));
end