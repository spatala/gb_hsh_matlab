function [Q] = sp_orth(A)
% sp_orth finds an orthogonal basis for the column space of A, but is 
% actually just a wrapper for qr.
% 
% Inputs:
%   A - matrix specifying the column space.
%
% Outputs:
%   Q - matrix whose columns form an orthonormal column space for A.
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    [Q, ~] = qr(A, 0);
end