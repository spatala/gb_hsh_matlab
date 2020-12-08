function [Q] = sp_orth(A)
% 
% Finds the **left** and **right** null spaces of the sparse matrix A.
% 
% - Inputs:
%   + A  : a sparse matrix with `m` rows and `n` columns.
%   + opt: `0` or `1` specifies that the left or right null 
%             space is desired.
%   + TOL: threshold for a pivot to be considered significant.
% 
% - Outputs:
%   + N: the left or right null space of the matrix A depending on the
%         value of opt. The sparsity is increased at the expense of the 
%         columns not being orthonormal.
% 
% - Notes:
%   + The algorithm is not well-established in the literature, 
%       but performs the decomposition A = L U Q, 
%     - where L and Q are invertible matrices, 
% 	  - U and A are the same size, and 
% 	  - U has the block structure U = [U_{11}, 0; 0, 0]. 
%   + If U is rank r, then since A Q^{-1} = L U the n - r columns 
%       of Q^{-1} corresponding to pivotless columns of U are a basis 
%       for the right null space. 
%   + Similarly, since L^{-1} A = U Q the m - r rows of L^{-1}
%       corresponding to pivotless rows of U are a basis for the 
%       left null space.
% 
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    [Q, ~] = qr(A, 0);
end