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