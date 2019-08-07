function [A] = clean_new(A, TOL)
% [A] = clean(S, TOL) - attempts to construct a column space for A with the
%   minimal number of nonzero entries. The matrix A can be complex.

% Forward sweep
for a = 1:size(A, 2)
    a
    r1 = find(abs(A(:, a)) > TOL, 1);
    A(:, a) = A(:, a) * conj(A(r1, a));
    A(:, a) = A(:, a) / sqrt(A(:, a)' * A(:, a));
    for b = (a + 1):size(A, 2)
        r2 = find(abs(A(:, b)) > TOL, 1);
        if r2 < r1
            swap = A(:, a);
            A(:, a) = A(:, b);
            A(:, b) = swap;
            
            r1 = find(abs(A(:, a)) > TOL, 1);
            A(:, a) = A(:, a) * conj(A(r1, a));
            A(:, a) = A(:, a) / sqrt(A(:, a)' * A(:, a));
        elseif r2 == r1
            A(:, b) = A(:, b) - (A(r1, b) / A(r1, a)) * A(:, a);
            A(:, b) = A(:, b) / sqrt(A(:, b)' * A(:, b));
        end
    end
end

% Backward sweep
for a = fliplr(1:size(A, 2))
    a
    A(:, a) = A(:, a) * conj(A(find(abs(A(:, a)) > TOL, 1), a));
    A(:, a) = A(:, a) / sqrt(A(:, a)' * A(:, a));
    for b = fliplr(1:(a - 1))
        A(:, b) = A(:, b) - (A(:, a)' * A(:, b)) * A(:, a);
        A(:, b) = A(:, b) / sqrt(A(:, b)' * A(:, b));
    end
end

% Attempt to remove remainders
Re = real(A);
Re(abs(Re) < TOL) = 0;
Im = imag(A);
Im(abs(Im) < TOL) = 0;
A = Re + 1i * Im;

% First nonzero coefficient of the eigenvector is postive.
for a = 1:size(A, 2)
    a
    A(:, a) = sign(A(find(abs(A(:, a)) > TOL, 1), a)) * A(:, a);
end
