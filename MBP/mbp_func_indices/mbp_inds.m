function tot_inds = mbp_inds(N)
tot_num = 0;
for ct1=0:N
    for ct2=0:N
        ct3 = min(ct1, ct2);
        tot_num = tot_num + (2*ct1+1)*(2*ct2+1)*(2*ct3+1);
    end
end
tot_inds = zeros(tot_num,5);
ind_start = 1;
for ct1=0:N
    for ct2=0:N
        inds_mat = mbp_ab(ct1,ct2);
        ind_stop = ind_start+size(inds_mat,1)-1;
        tot_inds(ind_start:ind_stop,1) = ind_start:ind_stop;
        tot_inds(ind_start:ind_stop,2:5) = inds_mat;
        ind_start = ind_stop+1;
    end
end
end