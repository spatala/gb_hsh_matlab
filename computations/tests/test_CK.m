% text_CK verifies that the conversion of a rotation angle and spherical angles
% of the rotation axis to a Cayley-Klein matrix and back returns the original
% result. This is not in general sufficient though, since a change in sign of
% an irrep of SU(2) (this occurs when multiplying rotations) can shift the
% apparent rotation angle outside of the allowed intervals.
%
% Requires access to the angles_to_CK.m and CK_to_angles.m files.
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.

% rng(142857);

w1  = pi * rand();
th1 = pi * rand();
ph1 = 2. * pi * rand();

CK = angles_to_CK(w1, th1, ph1);
[w2, th2, ph2] = CK_to_angles(CK);

disp(['Max difference: ', num2str([w1, th1, ph1] - [w2, th2, ph2])]);
