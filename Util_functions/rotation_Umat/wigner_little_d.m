function [d, err] = wigner_little_d(j, theta)
%
% Computes the Wigner d matrix for the given order and rotation angle, 
% with rows and colums labelled in increasing values of m.
% 
% - Input
%    + j: order of the matrix
%    + theta: rotation angle
% 
% 
% - Output
%    + d: Wigner (small) d-matrix
%    + err: estimate of the magnitude of numerical error.
% 
% - Notes  
%   + Follows the approach of X. M. Feng et al in 10.1103/PhysRevE.92.043307.
%   

m = -j:j;
n = 2 * j + 1;

X = sqrt((j + m) .* (j - m + 1)) / (2 * 1i);
Jy = zeros(n);
for a = 2:n
    b = n - a + 1;
    Jy(a - 1, a) = -X(a);
    Jy(b + 1, b) =  X(a);
end
% Jy is a tridiagonal Hermitian matrix
% Requires that eigenvectors be ordered with increasing eigenvalues
[V, ~] = eig(Jy);

W = V;
for a = 1:n
    W(:, a) = W(:, a) * exp(-1i * m(a) * theta);
end

d = W * V';

if nargout > 1
    err = max(abs(imag(d(:))));
end
d = real(d);
end