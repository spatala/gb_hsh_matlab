function [U] = rotation_mat(j, rot_angs)
%
% Input Parameters
%   j: Umatrix of order j.
%   rot_angs: 3 X 1 array.
%                w = rot_angs(1); th = rot_angs(2); ph = rot_angs(3);
%       + w: rotation angle;
%       + th: polar angle of rotation axis
%       + ph: azimuthal angle of rotation axis
%
% Returns
%   [U]: rotation matrix U of the order j. Size: (2*j+1, 2*j+1)
%
%   Follows Eq. 6 on page 81 of D. A. Varshalovich et al, Quantum Theory of
%   Angular Momentum, 1988.
%

w = rot_angs(1); th = rot_angs(2); ph = rot_angs(3);
if (abs(mod(w,2*pi)) < 1e-14)
    U = eye(2*j+1,2*j+1);
else
    tmp = tan(w / 2.) * cos(th);
    tmp = (1. - 1i * tmp) / realsqrt(1 + tmp^2.);
    r_base = 1i * exp(-1i * ph) * tmp;
    c_base = -1i * exp(1i * ph) * tmp;
    
    % Require $w \in [-\pi, \pi]$
    w = mod(w + pi, 2. * pi) - pi;
    xi = 2. * asin(sin(w / 2.) * sin(th));
    U = wigner_little_d(j, xi);
    
    m = -j:j;
    n = 2 * j + 1;
    for a = 1:n
        U(a, :) = U(a, :) * r_base^m(a);
        U(:, a) = U(:, a) * c_base^m(a);
    end
end
end