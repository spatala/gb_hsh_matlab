function tot_inds = mbp_inds(N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Get indices for MBP functions all the way to a=b=N.
%%%%%%%
%%%%%%% Input
%%%%%%% N: Integer
%%%%%%%         The maximum order for MBP functions.
%%%%%%%         We impose a = b = N.
%%%%%%%
%%%%%%% Output
%%%%%%% mbp: N X 7 matrix.
%%%%%%%         The columns are as follows:
%%%%%%%         1) Index number for fixed MBP function.
%%%%%%%         2) Index number for fixed a, b.
%%%%%%%         3) Value of 'a'
%%%%%%%         4) Value of 'b'
%%%%%%%         5) Value of 'gamma' - range [-c, c], where $c = min(a,b)$
%%%%%%%         6) Value of 'alpha' - range [-a, a]
%%%%%%%         7) Value of 'beta'  - range [-b, b]
%%%%%%%
tot_num = 0;
for ct1=0:N
    for ct2=0:N
        ct3 = min(ct1, ct2);
        tot_num = tot_num + (2*ct1+1)*(2*ct2+1)*(2*ct3+1);
    end
end
tot_inds = zeros(tot_num,7);
ind_start = 1;
for ct1=0:N
    for ct2=0:N
        inds_mat = mbp_ab(ct1,ct2);
        ind_stop = ind_start+size(inds_mat,1)-1;
        tot_inds(ind_start:ind_stop,1) = ind_start:ind_stop;
        tot_inds(ind_start:ind_stop,2:7) = inds_mat;
        ind_start = ind_stop+1;
    end
end
end