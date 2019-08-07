% text_R verifies that the irreps of SO(4) are unchanged by breaking apart the
% constituent rotations and then combining several of the irreps of SO(3) by
% the Clebsch-Gordan expansion. This is one of the essential steps in the 
% derivation of the simplified form of the basis functions in the MBP 
% parameterization. Equation numbers refer to particular equations of the
% manuscript, and will likely change in the future.
%
% Requires access to angles_to_CK.m, CK_to_angles.m, rotation_mat.m, and 
% clebsch_gordan.m.
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.

% clear;
% rng(142857);

w_m  = 8 * pi * (rand() - 1);
th_m = pi * rand();
ph_m = 2. * pi * rand();

w_b  = 8 * pi * (rand() - 1);
th_b = pi / 2;
ph_b = 2. * pi * rand();

w_z  = 8 * pi * (rand() - 1);
th_z = 0.;
ph_z = 0.;

CK_m = angles_to_CK(w_m, th_m, ph_m);
CK_b = angles_to_CK(w_b, th_b, ph_b);
CK_z = angles_to_CK(w_z, th_z, ph_z);

[w_1, th_1, ph_1] = CK_to_angles(CK_z * CK_b);
[w_2, th_2, ph_2] = CK_to_angles(CK_z * CK_b * CK_m);

a = 2;
b = 1;

% Equation 5
% Rows ordered by (\alpha', \beta') in lexicographic order
R1 = kron(rotation_mat(a, -w_1, th_1, ph_1), rotation_mat(b, -w_2, th_2, ph_2));
R1 = R1 * sqrt((2 * a + 1) * (2 * b + 1)) / (2 * pi^2);

% First equation on page 13
% Rows ordered by (\alpha', \beta') in lexicographic order
U1 = rotation_mat(a, -w_b, th_b, ph_b) * rotation_mat(a, -w_z, th_z, ph_z);
U2 = rotation_mat(b, -w_m, th_m, ph_m) * rotation_mat(b, -w_b, th_b, ph_b) * rotation_mat(b, -w_z, th_z, ph_z);
R2 = kron(U1, U2) * sqrt((2 * a + 1) * (2 * b + 1)) / (2 * pi^2);

disp(['Max R1 R2 difference: ', num2str(max(abs(R1(:) - R2(:))))]);

% Second equation on page 13
alpha = -a:a;
na = 2 * a + 1;

beta  = -b:b;
nb = 2 * b + 1;

% Rows ordered by (\alpha', \beta') in lexicographic order
R3 = zeros(na * nb, na * nb);

U_m = rotation_mat(b, -w_m, th_m, ph_m);
for e = abs(a - b):(a + b)
    epsilon = -e:e;
    ne = 2 * e + 1;
    U_e = rotation_mat(e, -w_b, th_b, ph_b);
    for p = 1:ne % \epsilon'
        [C, m1, m2] = clebsch_gordan(a, b, e, epsilon(p));
        ind = (m1 + a) * nb + m2 + b + 1;
        Cp = zeros(na * nb, 1);
        Cp(ind) = C;
        for q = 1:ne % \epsilon
            [C, m1, m2] = clebsch_gordan(a, b, e, epsilon(q));
            ind = (m1 + a) * nb + m2 + b + 1;
            Ce = zeros(na * nb, 1);
            Ce(ind) = C;
            for r = 1:nb % \delta'
                for s = 1:na % \alpha'
                    for t = 1:nb % \beta'
                        for u = 1:na % \alpha
                            for v = 1:nb % \beta
                                row = (alpha(s) + a) * nb + beta(t) + b + 1;
                                col = (alpha(u) + a) * nb + beta(v) + b + 1;
                                R3(row, col) = R3(row, col) + ...
                                    Cp((alpha(s) + a) * nb + beta(r) + b + 1) * ...
                                    Ce(col) * ...
                                    U_m(t, r) * ...
                                    U_e(p, q);
                            end
                        end
                    end
                end
            end
        end
    end
end

for p = 1:na
    for q = 1:nb
        col = (alpha(p) + a) * nb + beta(q) + b + 1;
        R3(:, col) = R3(:, col) * exp(1i * (alpha(p) + beta(q)) * w_z);
    end
end
R3 = R3 * sqrt((2 * a + 1) * (2 * b + 1)) / (2 * pi^2);

disp(['Max R1 R3 difference: ', num2str(max(abs(R1(:) - R3(:))))]);
