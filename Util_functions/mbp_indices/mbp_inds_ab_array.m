function tot_inds = mbp_inds_ab_array(symm_orders)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Get indices for MBP function with values 
%%%%% of (a,b) given by symm_orders.
%%%%%
%%%%% Input
%%%%% symm_orders: n X 2 array.
%%%%%         The possible values of (a,b) are provided in the symm_orders
%%%%%         array.
%%%%%
%%%%% Output
%%%%% mbp: N X 7 matrix.
%%%%%         The columns are as follows:
%%%%%         1) Index number for fixed MBP function.
%%%%%         2) Index number for fixed a, b.
%%%%%         3) Value of 'a'
%%%%%         4) Value of 'b'
%%%%%         5) Value of 'gamma' - range [-c, c], where $c = min(a,b)$
%%%%%         6) Value of 'alpha' - range [-a, a]
%%%%%         7) Value of 'beta'  - range [-b, b]
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a1 = symm_orders(:,1); b1 = symm_orders(:,2);  c1 = min(a1,b1);
num_rows = sum((2*a1+1).*(2*b1+1).*(2*c1+1));

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