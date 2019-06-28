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
    c_val = min(a_val, b_val);
    Na = 2*a_val; Nb = 2*b_val; Nc = 2*c_val; nsz = (Na+1)*(Nb+1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% Only works (so far) for
    %%%%% 1) X,Y,Z-pi rotations!; 2) Z-theta rotations!
    Ypi = [0,1,0,pi]; gl1 = vrrotvec2mat(Ypi); gl2 = vrrotvec2mat(Ypi);
    % Zth = [0,0,1,rand()]; gl1 = vrrotvec2mat(Zth); gl2 = vrrotvec2mat(Zth);
    
    gr1 = eye(3); gr2 = eye(3);
    Rr_ab_12 = so4_irrep(gr1,gr2,Na,Nb);
    Rl_ab_12 = so4_irrep(gl1,gl2,Na,Nb);
    trRr_ab_12 = transpose(Rr_ab_12);
    tr_Rrot = kron(Rl_ab_12,trRr_ab_12);
    Mrot = calc_Mrot_mat_Rab(tr_Rrot, a_val, b_val);
    
    Mrot(abs(Mrot)<1e-10) = 0; Mrot = sparse(Mrot);
    
    nsz1 = nsz*(Nc+1);
    ind_stop = ind_start + nsz1 - 1;
    %     [ind_start, ind_stop]
    Ypi_symm_mat(ind_start:ind_stop, ind_start:ind_stop) = Mrot;
    ind_start = ind_stop + 1;
end
end