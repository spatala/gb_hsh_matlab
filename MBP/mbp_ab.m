function inds_mat = mbp_ab(a,b)
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