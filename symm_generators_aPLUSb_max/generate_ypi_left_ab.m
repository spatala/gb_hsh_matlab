function ypileft_mat = generate_ypi_left_ab(a,b)
% 
% The matrix equivalent for (g1, g2) -> (Ypi*g1, Ypi*g2)
% M_ab(Ypi*g1, Ypi*g2) = M_ab(g1,g2) * ypileft_mat
% 
% Input:
% a,b:
%       Integers. Order for M_{a,b} function.
% 
% Output:
% ypileft_mat:
%       Matrix operation with size (Na+1)(Nb+1)(Nc+1)
%       Na = 2a; Nb = 2b; Nc = 2*min(a,b)
% 
c = min(a,b);
num_cols = (2*a+1)*(2*b+1)*(2*c+1);
tot_inds = mbp_ab(a,b);
tconst = (-1)^(a+b);
ypileft_mat = spalloc(num_cols, num_cols, num_cols);
for ct1=1:num_cols
    gamma1 = tot_inds(ct1,4);
    alpha1 = tot_inds(ct1,5);
    beta1 = tot_inds(ct1,6);
    gamma2 = -gamma1;
    %%%% Find the index for (-gamma,beta, alpha) row.
    ind1 = find(...
        (tot_inds(:,4) == gamma2 ) & ...
        (tot_inds(:,5) == alpha1 ) & ...
        (tot_inds(:,6) == beta1  ));
    %%%% The (a,b, gamma, alpha, beta)th column contains a (-1)^(a+b)
    %%%% in the (a,b, -gamma, beta, alpha)th row
    ypileft_mat(ind1,ct1) = tconst;
    %     ind_maps(ct1,:) = [ct1, ind1]; 
end
end

