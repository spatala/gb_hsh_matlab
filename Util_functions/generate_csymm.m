function csymm_mat = generate_csymm(gs1, gs2, symm_orders)

a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1,b1);
num_cols = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
nsymm = size(symm_orders,1);

csymm_mat = sparse(num_cols,num_cols);
ind_start = 1;

for ct1 = 1:nsymm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a_val = symm_orders(ct1,1); b_val = symm_orders(ct1,2);
    c_val = min(a_val, b_val);
    Na = 2*a_val; Nb = 2*b_val; Nc = 2*c_val; 
    nsz = (Na+1)*(Nb+1); nsz1 = nsz*(Nc+1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Rr_ab_12 = so4_irrep(gs1,gs2,Na,Nb); trRr_ab_12 = transpose(Rr_ab_12);
    Mrot = kron(eye(Nc+1), trRr_ab_12);
%     [a_val, b_val, c_val, size(Mrot,1)]
    
    ind_stop = ind_start + nsz1 - 1;
%     [ind_start, ind_stop]
    csymm_mat(ind_start:ind_stop, ind_start:ind_stop) = Mrot;
    ind_start = ind_stop + 1;
end

% mat_name = ['CS_Nmax_',num2str(Nmax),'.mat']; 
% save(mat_name,'csymm_mat');
end