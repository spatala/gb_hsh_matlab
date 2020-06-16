function [M] = mbp_basis(a, b, mbp_angs)
%
% Returns the basis functions $M^{a,b}$ for the grain boundary space
% using the MBP parameterization.
% 
% - Input:
%   + a, b: Indices for $M^{a,b}$ function.
%   + mbp_angs: (omega_m, theta_m, phi_m, omega_b, phi_b)
% 
% - Output:
%   + M: MBP basis function with order (a,b). 
% 	- Rows ordered by (\gamma, \alpha, \beta) in lexicographic order, i.e., starting with negative values, ending with positive values).
% 
% - Notes:
%   + Follows Equation () of the manuscript (finalize eqn. number after publishing).
%   + Not vectorized for array of mbp_angs.
%

t = num2cell(mbp_angs);
[w_m,th_m,ph_m,w_b,ph_b]=t{:};

CK_m = angles_to_CK(w_m, th_m, ph_m);
CK_b = angles_to_CK(w_b, pi / 2., ph_b);

[w_1, th_1, ph_1] = CK_to_angles(CK_b);
[w_2, th_2, ph_2] = CK_to_angles(CK_b * CK_m);

U1 = rotation_mat(a, [-w_1, th_1, ph_1]);
U2 = rotation_mat(b, [-w_2, th_2, ph_2]);

na = 2 * a + 1;
nb = 2 * b + 1;

gamma = -min(a, b):min(a, b);
ng = 2 * min(a, b) + 1;

M = zeros(na * nb, ng);
for p = 1:ng % \gamma
    M(:, p) = kron(U1(:, a + 1 + gamma(p)), U2(:, b + 1 - gamma(p)));
end
M = M * realsqrt((2 * a + 1) * (2 * b + 1)) / 7.874804972861210;

M = M(:);
end