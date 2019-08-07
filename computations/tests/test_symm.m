% test_symm verifies that the values of the symmetrized hyperspherical
% harmonics are the same (within numerical error) for symmetrically equivalent
% rotations. Also useful as an example of how to evaluate these functions.
% 
% Inputs:
%   N   - upper index of the hyperspherical harmonics.
%   generators - generators for the crystal point symmetry group as a cell
%         array of Cayley-Klein matrices. These should be identical to the
%         generators used when calculating the coefficients with make_symm.
%   trials - number of random rotations to test. Every one of these is right
%         multiplied by each element of the crystal point symmetry group, and
%         the resulting values are compared.
%   TOL - tolerance for a coefficient to be significant enough to include in
%         the definition of the symmetrized hyperspherical harmonics.
%
% Outputs:
%   check_diff - measures the spread in the function values for each of the
%         trial rotations. Includes variation with elements of the symmetry
%         group for all the symmetrized hyperspherical harmonics.
% 
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.

%#ok<*CLALL>
%#ok<*SAGROW>

clear all;
% rng(142857);

N = 8;

generators = {angles_to_CK(pi / 2., 0., 0.), ...
              angles_to_CK(pi / 2., pi / 2., 0.)};

trials = 100;

TOL = sqrt(eps);

% Generate symmetries
len0 = 0;
while true
    len1 = length(generators);
    for a = 1:len1
        for b = max(a, len0 + 1):len1
            CK = generators{a} * generators{b};
            include = true;
            for c = 1:length(generators)
                if norm(CK - generators{c}) < TOL || ...
                   norm(CK + generators{c}) < TOL
                    include = false;
                    break;
                end
            end
            if include
                generators{end + 1} = CK;
            end
        end
    end
    if length(generators) == len1
        break;
    else
        len0 = len1;
    end
end

% Generate test cases
w_t  = 8. * pi * rand(trials, 1);
th_t = 8. * pi * rand(trials, 1);
ph_t = 8. * pi * rand(trials, 1);

% Load coefficients
symm_coeff = load([num2str(N), '.mat']);
X = symm_coeff.X;
L = symm_coeff.L;
M = symm_coeff.M;

% Calculate values and check spread
check_diff = zeros(trials, 1);
f_vals = zeros(numel(generators), size(X, 2));
for a = 1:trials
    CK = angles_to_CK(w_t(a), th_t(a), ph_t(a));
    for b = 1:numel(generators)
        [w, th, ph] = CK_to_angles(CK * generators{b});
        for c = 1:size(X, 2)
            for d = 1:(N + 1)^2
                if abs(X(d, c)) > TOL
                    f_vals(b, c) = f_vals(b, c) + X(d, c) * ...
                        get_hsh(N, L(d), M(d), w, th, ph);
                end
            end
        end
    end
    check_diff(a) = norm(max(f_vals) - min(f_vals));
end

disp(['Max error: ', num2str(max(check_diff))]);
