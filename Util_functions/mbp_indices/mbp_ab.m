function inds_mat = mbp_ab(a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Get indices for MBP functions with fixed a and b.
%%%%%%%
%%%%%%% Input
%%%%%%% a, b: Integers
%%%%%%%         order of MBP function (gamma_MBP^{a,b}_{alpha, beta})
%%%%%%%
%%%%%%% Output
%%%%%%% mbp: N X 6 matrix.
%%%%%%%         The columns are as follows:
%%%%%%%         1) Index number for fixed a, b.
%%%%%%%         2) Value of 'a'
%%%%%%%         3) Value of 'b'
%%%%%%%         4) Value of 'gamma' - range [-c, c], where $c = min(a,b)$
%%%%%%%         5) Value of 'alpha' - range [-a, a]
%%%%%%%         6) Value of 'beta'  - range [-b, b]
%%%%%%%

c = min(a,b);
tot_num = (2*a+1)*(2*b+1)*(2*c+1);
inds1 = 1:tot_num;
ct1 = 1;
inds_mat = zeros(tot_num,6);
for g1=-c:c
    for a1=-a:a
        for b1=-b:b
            inds_mat(ct1, :) = [inds1(ct1), a, b, g1, a1, b1];
            ct1=ct1+1;
        end
    end
end

end