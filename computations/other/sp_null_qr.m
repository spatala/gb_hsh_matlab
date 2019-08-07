function [null_space] = sp_null_qr(A, opt, TOL)
% sp_null_qr finds the left and right null spaces of the sparse matrix A. The 
% algorithm is not well-established in the literature, but performs the
% decomposition A = Q R E where Q and E are invertible matrices, R and A are
% the same size, and R has the block structure R = [R_{11}, 0; 0, 0]. If R is
% rank r, then since A E^{-1} = Q R the n - r columns of E^{-1} corresponding
% to pivotless columns of R are a basis for the right null space. Similarly,
% since Q^{-1} A = R E the m - r rows of Q^{-1} corresponding to pivotless
% rows of R are a basis for the left null space.
% 
% Inputs:
%   A   - a sparse matrix with m rows and n columns.
%   opt - 0 or 1 specifies that the left or right null space is desired.
%   TOL - threshold for a pivot to be considered significant.
%
% Outputs:
%   null_space - the left or right null space of the matrix A depending on the
%         value of opt. The sparsity is increased at the expense of the columns
%         not being orthonormal.
% 
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    switch opt
        case 0
            I = find(A);
            A(I(abs(A(I)) < TOL)) = 0.;
            [Q, R, ~] = qre(A, TOL);
            r = nnz(abs(spdiags(R, 0)) > TOL);
            inv_Q = inv(Q);
            null_space = inv_Q((r + 1):size(A, 1), :);
        case 1
            null_space = sp_null2(A.', 0, TOL).';
        otherwise
            disp('Invalid option. Consult documentation.');
    end
end

function [Q, R, E] = qre(A, TOL)
% luq is the workhorse of sp_null2, and actually performs the decomposition
% A = Q * R * E. This is constructed using the A = Q * R decomposition
% given by the qr function in MATLAB.
% 
% Inputs:
%   A   - a sparse matrix with m rows and n columns.
%   TOL - threshold for a pivot to be considered significant.
%
% Outputs:
%   Q   - an invertible matrix with the same number of rows as A.
%   R   - a matrix with the block structure R = [R_{11}, 0; 0, 0] where 
%         R_{11} is a square matrix indicating the rank of A.
%   E   - an invertible matrix with the same number of columns as A.
% 
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    [m, n] = size(A);
    if ~issparse(A)
        A = sparse(A);
    end    

    % This should be as stable as possible
    [Q, R] = qr(A);
    E = speye(n);
    
    % Find the row and col indices of significant pivots
    p1 = find(abs(spdiags(R, 0)) > TOL)';
    p = numel(p1);
    
    % Permute rows and columns to give R = [R_{11}, R_{12}; R_{21}, R_{22}]
    % with all of the significant pivots in R_{11}.
    r2 = 1:m;
    r2(p1) = [];
    r_perm = [p1, r2];
    Q = Q(:, r_perm);
    R = R(r_perm, :);
    
    c2 = 1:n;
    c2(p1) = [];
    c_perm = [p1, c2];
    R = R(:, c_perm);
    E = E(c_perm, :);
    
    % Update indices to reflect reordering
    p1 = 1:p;
    r2 = (p + 1):m;
    c2 = (p + 1):n;
    
    % The R_{21} block is converted to zeros by R = [I, 0; -X, I] * R and 
    % Q = Q * [I, 0; X, I], preserving the product A = Q * R * E.
    X = R(r2, p1) / R(p1, p1);
    Q(:, p1) = Q(:, p1) + Q(:, r2) * X;
    R(r2, p1) = 0.;
    R(r2, c2) = R(r2, c2) - X * R(p1, c2);
    
    % The R_{12} block is converted to zeros by R = R * [I, -X; 0, I] and 
    % E = [I, X; 0, I] * E, preserving the product A = Q * R * E. R_{22} is
    % unchanged as a consequence of the previous step.
    X = R(p1, p1) \ R(p1, c2);
    R(p1, c2) = 0.;
    E(p1, :) = E(p1, :) + X * E(c2, :);
    
    % At this point R = [R_{11}, 0; 0, R_{22}] but R_{22} could still have
    % nonzero values off of the diagonal. Permute rows and columns to
    % change block structure to R = [R_{11}, 0, 0; 0, R_{22}, 0; 0, 0, 0].
    r3 = p + find(max(abs(R(r2, c2)), [], 2) > TOL)';
    c3 = p + find(max(abs(R(r2, c2)), [], 1) > TOL);
    if ~isempty(r3) || ~isempty(c3)
        r4 = r2;
        r4(r3 - p) = [];
        r_perm = [r3, r4];
        Q(:, (p + 1):m) = Q(:, r_perm);
        R((p + 1):m, :) = R(r_perm, :);

        c4 = c2;
        c4(c3 - p) = [];
        c_perm = [c3, c4];
        R(:, (p + 1):n) = R(:, c_perm);
        E((p + 1):n, :) = E(c_perm, :);
        
        % Update indices to reflect reordering
        r3 = p + (1:numel(r3));
        c3 = p + (1:numel(c3));
        
        % Perform decomposition of R_{22} = Q_r * R_r * E_r, resulting in
        % A = [Q_1, Q_2 * Q_r, Q_3] * ...
        %     [R_{11}, 0, 0; 0, R_r, 0; 0, 0, 0] * ...
        %     [E_1; E_r * E_2; E_3]
        [Qr, Rr, Er] = qre(R(r3, c3), TOL);
        R(r3, c3) = Rr;
        Q(:, r3) = Q(:, r3) * Qr;
        E(c3, :) = Er * E(c3, :);
        
        % Find the number of significant pivots
        q = nnz(abs(spdiags(R, 0)) > TOL);
        R = [R(1:q, 1:q), sparse(q, n - q); sparse(m - q, q), sparse(m - q, n - q)];
    else
        R(r2, c2) = 0.;
    end
end
