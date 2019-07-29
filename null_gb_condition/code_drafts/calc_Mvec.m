function Mvec = calc_Mvec(g1,g2,symm_orders)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Compute Mvec full vector given by symm_orders
%%%%%% 
%%%%%% Inputs: g1, g2, symm_orders
%%%%%% 
%%%%%% Output: Mvec
%%%%%% 
%%%%%% TODO: Have to make this code efficient
%%%%%%          1) Vectorize?
%%%%%%          2) No need to compute full R_ab
%%%%%%          

a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1,b1);
num_cols = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
nsymm = size(symm_orders,1);

Mvec = zeros(1,num_cols);
ind_start = 1;
for ct1=1:nsymm
    a_val = symm_orders(ct1,1); b_val = symm_orders(ct1,2);
    c_val = min(a_val, b_val);
    Na = 2*a_val; Nb = 2*b_val; Nc = 2*c_val; nsz = (Na+1)*(Nb+1);
    
    R_ab_12 = so4_irrep(g1,g2,Na,Nb);
    Mvec_ab_12 = calc_Mfunc_Rab(R_ab_12, a_val, b_val);
    
    nsz1 = nsz*(Nc+1);
    ind_stop = ind_start + nsz1 - 1;
    Mvec(ind_start:ind_stop) = Mvec_ab_12;
    ind_start = ind_stop + 1;
end

end