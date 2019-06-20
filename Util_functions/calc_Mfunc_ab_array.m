function M_vals = calc_Mfunc_ab_array(symm_orders, num_rows, rots)
nsymm = size(symm_orders,1);

tot_inds = zeros(num_rows,7);
M_vals = zeros(1,num_rows);
ind_start = 1;
for ct1=1:nsymm
    a_val = symm_orders(ct1,1);
    b_val = symm_orders(ct1,2);
    M = calc_Mfunc(a_val,b_val,rots);
    inds_mat = mbp_ab(a_val,b_val);
    ind_stop = ind_start+size(inds_mat,1)-1;
    tot_inds(ind_start:ind_stop,1) = ind_start:ind_stop;
    tot_inds(ind_start:ind_stop,2:7) = inds_mat;
    M_vals(ind_start:ind_stop) = M';
    ind_start = ind_stop+1;
end

end