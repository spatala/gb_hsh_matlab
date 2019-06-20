function tot_inds = mbp_inds_ab_array(symm_orders, num_rows)
tot_inds = zeros(num_rows,7);
ind_start = 1;

for ct1=1:size(symm_orders,1)
    a_val = symm_orders(ct1,1);
    b_val = symm_orders(ct1,2);
    
    inds_mat = mbp_ab(a_val,b_val);
    ind_stop = ind_start+size(inds_mat,1)-1;
    tot_inds(ind_start:ind_stop,1) = ind_start:ind_stop;
    tot_inds(ind_start:ind_stop,2:7) = inds_mat;
    ind_start = ind_stop+1;
end
end