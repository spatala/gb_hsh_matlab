function [M] = mbp_basis(a, b, w_m, th_m, ph_m, w_b, ph_b)
% mbp_basis calculates the basis functions on the grain boundary space using
% the misorientation/boundary plane convention.
% 
% Inputs:
%   a    - left upper index.
%   b    - right upper index.
%   w_m  - misorientation rotation angle.
%   th_m - polar angle of misorientation rotation axis.
%   ph_m - azimuthal angle of misorientation rotation axis.
%   w_m  - boundary plane rotation angle.
%   ph_m - azimuthal angle of boundary plane rotation axis.
%
% Outputs:
%   M    - column vector of all of the basis functions evaluated for the 
%          specified grain boundary. Rows labelled by (\gamma, \alpha, \beta)
%          in lexicographic order (start with negative values).
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    CK_m = angles_to_CK(w_m, th_m, ph_m);
    CK_b = angles_to_CK(w_b, pi / 2., ph_b);
    
    [w_2, th_2, ph_2] = CK_to_angles(CK_b * CK_m);

    U1 = rotation_mat(a, -w_b, pi / 2., ph_b);
    U2 = rotation_mat(b, -w_2, th_2, ph_2);
    
    na = 2 * a + 1;
    nb = 2 * b + 1;
    
    gamma = -min(a, b):min(a, b);
    ng = 2 * min(a, b) + 1;
    
    % Uses the product of irreps of SO(3)
    M = zeros(na * nb, ng);
    for p = 1:ng % \gamma
        M(:, p) = kron(U1(:, a + 1 + gamma(p)), U2(:, b + 1 - gamma(p)));
    end

    % Numerical value equal to \sqrt{2 \pi^3}
    M = M(:) * realsqrt((2 * a + 1) * (2 * b + 1)) / 7.874804972861210;
end
