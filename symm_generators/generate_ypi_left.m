function Ypi_symm_mat = generate_ypi_left(symm_orders)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1,b1);
num_cols = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
nsymm_ab = size(symm_orders,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ypi_symm_mat = sparse(num_cols, num_cols);

ind_start = 1;
for ct1 = 1:nsymm_ab
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a_val = symm_orders(ct1,1); b_val = symm_orders(ct1,2);
    Mrot = generate_ypi_left_ab(a_val,b_val);
    c_val = min(a_val, b_val);
    Na = 2*a_val; Nb = 2*b_val; Nc = 2*c_val; nsz = (Na+1)*(Nb+1);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%% Only works (so far) for
%     %%%%% 1) X,Y,Z-pi rotations!; 2) Z-theta rotations!
%     Ypi = [0,1,0,pi]; gl1 = vrrotvec2mat(Ypi); gl2 = vrrotvec2mat(Ypi);
%     % Zth = [0,0,1,rand()]; gl1 = vrrotvec2mat(Zth); gl2 = vrrotvec2mat(Zth);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     gr1 = eye(3); gr2 = eye(3);
%     Rr_ab_12 = so4_irrep(gr1,gr2,Na,Nb);
%     Rl_ab_12 = so4_irrep(gl1,gl2,Na,Nb);
%     trRr_ab_12 = transpose(Rr_ab_12);
%     tr_Rrot = kron(Rl_ab_12,trRr_ab_12);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Mrot = calc_Mrot_mat_Rab(tr_Rrot, a_val, b_val);
%     Mrot(abs(Mrot)<1e-10) = 0; Mrot = sparse(Mrot);
    
    
%     norm(full(Mrot)-full(Mrot1))
    
    nsz1 = nsz*(Nc+1);
    ind_stop = ind_start + nsz1 - 1;
    Ypi_symm_mat(ind_start:ind_stop, ind_start:ind_stop) = Mrot;
    ind_start = ind_stop + 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end

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