function [D] = wigner_d(j, al, be, ga)
% wigner_d converts a triple of Euler angles (rotations about z, y, and z) into
% the equivalent irrep of SO(3). Beware that the precise interpretation of the
% Euler angles is not obvious, but converting the Euler angles into a rotation
% angle and spherical angles of the rotation axis using Eq. 16 on page 26 of 
% D. A. Varshalovich et al, Quantum Theory of Angular Momentum, 1988 and using
% rotation_mat gives an identical irrep of SO(3).
% 
% Inputs:
%   j  - specifies dimension of the irrep.
%   al - initial z rotation in the interval [0, 2 \pi].
%   be - initial y rotation in the interval [0, \pi].
%   ga - initial z rotation in the interval [0, 2 \pi].
%
% Outputs:
%   D  - (2 j + 1)-dimensional representation of SO(3), using the conventions
%        established in Eq. 1 on page 76 of D. A. Varshalovich et al, Quantum
%        Theory of Angular Momentum, 1988. Rows and columns ordered in
%        increasing values of m' and m.
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    D = wigner_little_d(j, be);

    m = -j:j;
    n = 2 * j + 1;
    for a = 1:n
        D(a, :) = D(a, :) * exp(-1i * m(a) * al);
        D(:, a) = D(:, a) * exp(-1i * m(a) * ga);
    end
end
