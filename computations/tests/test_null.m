% test_null verifies that three versions of the basis functions on the grain
% boundary space restricted to the null boundary (identity misorientation) are
% equivalent. The first uses mbp_basis and is intended as the general purpose
% function. The second uses a simplification of the second version of the irrep
% of SO(4) considered in test_R.m. The third uses the equation that appears in
% the section on the null boundary, and is conceptually identical to the second
% with the restricted matrix elements of the irrep of SO(3) replaced by
% spherical harmonics. Equation numbers refer to particular equations of the
% manuscript, and will likely change in the future.
%
% Requires access to mbp_basis.m, rotation_mat.m, clebsch_gordan.m, and
% spherical_harmonics.m.
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.

% clear;
% rng(142857);

w_m  = 0.;
th_m = pi * rand();
ph_m = 2. * pi * rand();

w_b  = 8 * pi * (rand() - 1);
th_b = pi / 2;
ph_b = 2. * pi * rand();

a = 6;
b = 3;

% Equation 14
% Rows ordered by (\gamma, \alpha, \beta) in lexicographic order
M1 = mbp_basis(a, b, w_m, th_m, ph_m, w_b, ph_b);

% Restriction of equation on page 13
% M^{a b}_{\alpha \beta \gamma}(g_b, g_b) = \sum_{e \epsilon} \frac{\Pi_{a
% b}}{\sqrt{2 \pi^3}} C^{e \epsilon}_{a \alpha b \beta} C^{e 0}_{a \gamma b
% -\gamma} U^{e}_{\epsilon 0}(g_b^{-1})
alpha = -a:a;
na = 2 * a + 1;

beta  = -b:b;
nb = 2 * b + 1;

gamma = -min(a, b):min(a, b);
ng = 2 * min(a, b) + 1;

% Rows ordered by (\alpha, \beta) in lexicographic order
M2 = zeros(na * nb, ng);
for e = abs(a - b):(a + b)
    epsilon = -e:e;
    ne = 2 * e + 1;
    U_e = rotation_mat(e, -w_b, th_b, ph_b);
    
    [C, m1, m2] = clebsch_gordan(a, b, e, 0);
    ind = (m1 + a) * nb + m2 + b + 1;
    Cz = zeros(na * nb, 1);
    Cz(ind) = C;
    for p = 1:ne % \epsilon
        [C, m1, m2] = clebsch_gordan(a, b, e, epsilon(p));
        ind = (m1 + a) * nb + m2 + b + 1;
        Ce = zeros(na * nb, 1);
        Ce(ind) = C;
        for s = 1:na % \alpha
            for t = 1:nb % \beta
                for u = 1:ng % \gamma
                    row = (alpha(s) + a) * nb + beta(t) + b + 1;
                    M2(row, u) = M2(row, u) + ...
                        Ce(row) * ...
                        Cz((gamma(u) + a) * nb - gamma(u) + b + 1) * ...
                        U_e(p, e + 1);
                end
            end
        end
    end
end
M2 = M2 * sqrt((2 * a + 1) * (2 * b + 1)) / realsqrt(2 * pi^3);

disp(['Max M1 M2 difference: ', num2str(max(abs(M1(:) - M2(:))))]);

% First equation on page 20
% Rows ordered by (\alpha, \beta) in lexicographic order
M3 = zeros(na * nb, ng);
for e = abs(a - b):(a + b)
    epsilon = -e:e;
    ne = 2 * e + 1;
    Y_e = spherical_harmonics(e, w_b, ph_b);
    
    [C, m1, m2] = clebsch_gordan(a, b, e, 0);
    ind = (m1 + a) * nb + m2 + b + 1;
    Cz = zeros(na * nb, 1);
    Cz(ind) = C;    
    for p = 1:ne % \epsilon
        [C, m1, m2] = clebsch_gordan(a, b, e, -epsilon(p));
        ind = (m1 + a) * nb + m2 + b + 1;
        Ce = zeros(na * nb, 1);
        Ce(ind) = C;
        for s = 1:na % \alpha
            for t = 1:nb % \beta
                for u = 1:ng % \gamma
                    row = (alpha(s) + a) * nb + beta(t) + b + 1;
                    M3(row, u) = M3(row, u) + ...
                        (-1i)^epsilon(p) / sqrt(2 * e + 1)* ...
                        Ce(row) * ...
                        Cz((gamma(u) + a) * nb - gamma(u) + b + 1) * ...
                        Y_e(p);
                end
            end
        end
    end
end
M3 = M3 * sqrt(2) * sqrt((2 * a + 1) * (2 * b + 1)) / pi;

disp(['Max M1 M3 difference: ', num2str(max(abs(M1(:) - M3(:))))]);
