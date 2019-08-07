% test_M verifies that the calculation of the basis functions on the grain 
% boundary space by mbp_basis is equivalent to a restriction of the elements of
% the irreps of SO(4). Equation numbers refer to particular equations of the
% manuscript, and will likely change in the future.
%
% Requires access to angles_to_CK.m, CK_to_angles.m, mbp_basis.m, and 
% rotation_mat.m.
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.

clear;
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
b = 4;

% Equation 14
% Rows ordered by (\gamma, \alpha, \beta) in lexicographic order
M1 = mbp_basis(a, b, w_m, th_m, ph_m, w_b, ph_b);

% Equation 5
% Rows ordered by (\alpha', \beta') in lexicographic order
R1 = kron(rotation_mat(a, -w_1, th_1, ph_1), rotation_mat(b, -w_2, th_2, ph_2));
R1 = R1 * sqrt((2 * a + 1) * (2 * b + 1)) / (2 * pi^2);

na = 2 * a + 1;
nb = 2 * b + 1;

gamma = -min(a, b):min(a, b);
g_ind = (gamma + a) * nb - gamma + b + 1;

M2 = R1(:, g_ind);
M2 = M2(:) * realsqrt(2. * pi);

disp(['Max M1 M2 difference: ', num2str(max(abs(M1(:) - M2(:))))]);
