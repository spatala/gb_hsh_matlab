function [] = generate_gb_const()
a = 5;
b = 3;
mat_ab = null_mat_ab(a,b);

end
function mat_ab = null_mat_ab(a,b)
c = min(a,b);
na = 2*a+1; 
nb = 2*b+1;
ng = 2*c+1;

e_range = abs(a-b):a+b; 
num_e = numel(e_range);
mat_ab = sparse(sum(2*e_range+1),na*nb*ng);

r_ct = 1;
for e_ct = 1:num_e
    e = e_range(e_ct);
    [Cz, ~, ~] = clebsch_gordan(a, b, e, 0);
    eps_range = -e:e;
    ne = 2 * e + 1;
    for eps_ct = 1:ne
        eps1 = eps_range(eps_ct);
        [Ce, m1, m2] = clebsch_gordan(a, b, e, -eps1);
        ind = (m1 + a) * nb + m2 + b + 1;
        C_val = Ce*(Cz.');
        C_val = C_val(:);
        ind1 = ind + (0:ng-1)*na*nb;
        ind1 = ind1(:);
        mat_ab(r_ct, ind1) = C_val/(sqrt(2*e+1));
        r_ct = r_ct + 1;
    end
end
mat_ab = mat_ab*sqrt(2)*sqrt((2*a+1)*(2*b+1))/pi;
end
