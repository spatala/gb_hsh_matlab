function ypileft_mat = generate_ypi_left_ab(a,b)
c = min(a,b);
num_cols = (2*a+1)*(2*b+1)*(2*c+1);
tot_inds = mbp_ab(a,b);
tconst = (-1)^(a+b);
ypileft_mat = sparse(num_cols, num_cols);

for ct1=1:num_cols
    a1 = tot_inds(ct1,2); b1 = tot_inds(ct1,3);
    gamma1 = tot_inds(ct1,4);
    alpha1 = tot_inds(ct1,5); beta1 = tot_inds(ct1,6);
    gamma2 = -gamma1;
    
    %%%% Find the index for (a,b,-gamma,beta, alpha) row.
    ind1 = find(...
        (tot_inds(:,2) == a1     ) & ...
        (tot_inds(:,3) == b1     ) & ...
        (tot_inds(:,4) == gamma2 ) & ...
        (tot_inds(:,5) == alpha1 ) & ...
        (tot_inds(:,6) == beta1  ));
    
    %%%% The (a,b, gamma, alpha, beta)th column contains a (-1)^(a+b)
    %%%% in the (a,b, -gamma, beta, alpha)th row
    ypileft_mat(ind1,ct1) = tconst;
    %     ind_maps(ct1,:) = [ct1, ind1];
    
end
end