% rng(142857);

ERR = 1e-10;
TOL = 1e-8;

% A = [1, 4; 2, 8; 5, 7; 1, 4];
% A = [A, A + rand(size(A)) * ERR];
% A = sparse(A);

A = sparse(4,4);

A = A + rand(size(A))*ERR;

[~, SpRight] = spspaces(A, 2, TOL);
null_space = clean(SpRight{1}(:, SpRight{3}), TOL);

err = abs(A * null_space);
disp(['Max error: ', num2str(max(err(:)))]);
disp(full(null_space));