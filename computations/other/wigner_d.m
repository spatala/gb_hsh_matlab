function [U] = wigner_d(j, al, be, ga)
    U = wigner_little_d(j, be);
    m = -j:j;
    n = 2 * j + 1;
    for a = 1:n
        U(a, :) = U(a, :) * exp(-1i * m(a) * al);
        U(:, a) = U(:, a) * exp(-1i * m(a) * ga);
    end
end