function flip_mat = generate_flip_mat(symm_orders)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1,b1);
num_cols = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% The variable tot_inds contains the mapping between the column-index
%%%% for the M-function and the (a,b,gamma,alpha, beta) values.
tot_inds = mbp_inds_ab_array(symm_orders, num_cols);

%%%% flip_mat is the $A$ matrix that specifies the grain-exchange-symmetry
%%%% for the coefficients.
flip_mat = sparse(num_cols, num_cols);
% ind_maps = zeros(num_cols,2);
for ct1=1:num_cols
    %%%% Get the (a,b,gamma,alpha, beta) corresponding to the column.
    a1 = tot_inds(ct1,3); b1 = tot_inds(ct1,4);
    gamma1 = tot_inds(ct1,5); 
    alpha1 = tot_inds(ct1,6); beta1 = tot_inds(ct1,7);
    %%%% Find the index for (b,a,-gamma,beta, alpha) row.
    a2 = b1; b2 = a1; gamma2 = -gamma1;
    alpha2 = beta1; beta2 = alpha1;
    ind1 = find(...
        (tot_inds(:,3) == a2     ) & ...
        (tot_inds(:,4) == b2     ) & ...
        (tot_inds(:,5) == gamma2 ) & ...
        (tot_inds(:,6) == alpha2 ) & ...
        (tot_inds(:,7) == beta2  ));
    
    %%%% The (a,b, gamma, alpha, beta)th column contains a 1 in the 
    %%%% (b,a, -gamma, beta, alpha)th row
    flip_mat(ind1,ct1) = 1;
%     ind_maps(ct1,:) = [ct1, ind1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s1 = load(['flip_inds_map_',num2str(Nmax),'.mat']);
% ind_map = s1.ind_map;
% max(max(abs(ind_maps - ind_map)))
% %%% Save the flip-matrix in .mat file.
% mat_name = ['flip_Nmax_',num2str(Nmax),'.mat']; 
% save(mat_name, 'flip_mat','ind_maps');

end