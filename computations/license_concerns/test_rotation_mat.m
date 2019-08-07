% text_rotation_mat verifies that the calculation of the irreps of SO(3) in
% rotation_mat.m is equivalent to the prior implementation in rotation.m (with
% the rows and columns reordered). The current version is considerably faster
% and more stable because there are no explicit factorials appearing in the
% formulas. THE PREVIOUS VERSION USES A DIFFERENT LICENSE, and to avoid
% license complications should perhaps not be included in a distribution.
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

%#ok<*FLUDLR>

% rng(142857);

w  = 8 * pi * (rand() - 0.5);
th = pi * rand();
ph = 2 * pi * rand();

j = 5;

% Explicit form for j = 1 following page 127 of D. A. Varshalovich et al,
% Quantum Theory of Angular Momentum, 1988.
% 
% cw = cos(w / 2.);
% sw = sin(w / 2.);
% ct = cos(th);
% st = sin(th);
% ref = [(cw - 1i * sw * ct)^2, ...
%        -1i * sqrt(2) * sw * st * exp(-1i * ph) * (cw - 1i * sw * ct), ...
%        -(sw * st * exp(-1i * ph))^2;
%        -1i * sqrt(2) * sw * st * exp(1i * ph) * (cw - 1i * sw * ct), ...
%        1 - 2 * sw^2 * st^2, ...
%        -1i * sqrt(2) * sw * st * exp(-1i * ph) * (cw + 1i * sw * ct);
%        -(sw * st * exp(1i * ph))^2, ...
%        -1i * sqrt(2) * sw * st * exp(1i * ph) * (cw + 1i * sw * ct), ...
%        (cw + 1i * sw * ct)^2];
% ref = flipud(fliplr(ref));

U = rotation_mat(j, w, th, ph);

R = rotation([cos(ph) * sin(th), sin(ph) * sin(th), cos(th)], w, 2 * j);
R = flipud(fliplr(R));

disp(['Max U R difference: ', num2str(max(abs(U(:) - R(:))))]);

tic;
for a = 1:1000
    U = rotation_mat(j, w, th, ph);
end
toc;

ax = [cos(ph) * sin(th), sin(ph) * sin(th), cos(th)];
tic;
for a = 1:1000
    R = rotation(ax, w, 2 * j);
end
toc;
