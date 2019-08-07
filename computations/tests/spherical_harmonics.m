function [Ylm] = spherical_harmonics(l, th, ph)
% spherical_harmonics constructs the spherical harmonics for the specified
% total angular momentum and every allowed z component of angular momentum.
% Follows the Condon-Shortley convention.
% 
% Inputs:
%   l   - total angular momentum.
%   th  - polar angle in the interval [0, \pi].
%   ph  - aximuthal angle in the interval [0, 2 \pi].
%
% Outputs:
%   Ylm - column vector of the values of the spherical harmonics for all
%         allowed values of the z component of angular momentum. Rows are
%         ordered in increasing values of m. Follows the convention established
%         in Eq. 1 on page 133 of D. A. Varshalovich et al, Quantum Theory of
%         Angular Momentum, 1988.
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    
    % Need to be careful about inputs
    th2 = mod(th, pi);
    ph2 = mod(ph, 2. * pi);
    
    % Construct associated Legendre functions of all orders
    Plm = legendre(l, cos(th2));
    for m = 0:l
        Plm(m + 1) = Plm(m + 1) / realsqrt(prod((l - m + 1):(l + m)));
    end

    % Construct spherical harmonics
    Ylm = [flipud(Plm(2:(l + 1))); Plm];
    for m = -l:-1
        Ylm(l + m + 1) = Ylm(l + m + 1) * (-1)^abs(m);
    end
    for m = -l:l
        Ylm(l + m + 1) = Ylm(l + m + 1) * exp(1i * m * ph2);
    end
    % Numerical constant equal to \sqrt{4 \pi}
    Ylm = Ylm * realsqrt(2 * l + 1) / 3.544907701811032;
    
    % The following phase adjustments derive from Eqs. 10 and 11 on page 141 of
    % D. A. Varshalovich et al, Quantum Theory of Angular Momentum, 1988.

    % Phase change when th outside allowed bounds
    n = round(abs(th - th2) / pi);
    if mod(n, 2) == 1
        Ylm = (-1)^l * Ylm;
    end
    
    % Phase change when ph outside allowed bounds
    n = round(abs(ph - ph2) / pi);
    if mod(n, 2) == 1
        Ylm = (-1)^(-l:l) .* Ylm;
    end
end
