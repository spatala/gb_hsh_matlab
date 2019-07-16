function basis_eq = check_equi_basis(X1, X2)
nsz = size(X2,2);
if (nsz == size(X1,2))
    diff_vecs = zeros(nsz,1);
    tA1 = X1;
    for ct1 = 1:nsz
        tB1 = X2(:,ct1);
        diff_vecs(ct1) = norm(tA1*(tA1\tB1) - tB1);
    end
    if (norm(diff_vecs)/nsz < 1e-14)
       basis_eq = 1; 
    else
        basis_eq = 0;
    end
else
    basis_eq = 0;
end

end