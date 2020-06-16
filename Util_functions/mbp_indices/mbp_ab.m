function inds_mat = mbp_ab(a,b)
%
% Get indices for MBP functions with fixed a and b.
% 
% - Input
%   + a, b: Integers
%         order of MBP function 
%         ($\tensor[_{\gamma}]{M}{^a_{\alpha}^b_{\beta}}$)
% 
% - Output
%   + mbp: N X 6 matrix. The columns are as follows:
%     - Index number for fixed a, b.
%     - Value of 'a'
%     - Value of 'b'
%     - Value of 'gamma' - range [-c, c], where $c = min(a,b)$
%     - Value of 'alpha' - range [-a, a]
%     - Value of 'beta'  - range [-b, b]
%

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