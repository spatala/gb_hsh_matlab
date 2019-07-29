function Z = sph_harm_Mab_formula(ab_ord, l, th, phi)

Yl = compute_complex_Yl(th,phi,l);
m1 = compute_mult_factor(ab_ord, l);
Z = m1*Yl/sqrt(2*l+1);
end

function mult1 = compute_mult_factor(ab_ord,l)

a1 = ab_ord(1); b1 = ab_ord(2); g1 = ab_ord(3); al1 = ab_ord(4); 
be1 = ab_ord(5);

mult1 = zeros(1,2*l+1);
ct1 = 1;
for m=-l:l
    cleb1 = clebschgordan(a1, al1, b1, be1, l, -m);
    cleb2 = clebschgordan(a1, g1, b1, -g1, l, 0);
    %     [cleb1, cleb2]
    mult1(ct1) = ((-1i)^m)*cleb1*cleb2;
    ct1 = ct1 + 1;
end

end