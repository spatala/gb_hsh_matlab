% This script takes about 10 seconds to run. It is worthwhile changing the
% tolerances and seeing how the error and number of elements close to zero
% change.

TOL1 = 1e-14;
TOL2 = 1e-12;

load('null_mat_Nmax_10.mat')

% Calculates the full right null space
tic;
ns = sp_null(null_mat, 1, TOL1);
T = toc;
err = max(max(abs(null_mat * ns)));
disp(['Found null space in ', num2str(T), ...
      ' seconds with error of ', num2str(err)]);
disp(['Fractional fill is ', num2str(nnz(ns) / numel(ns))]);
disp(['Number of elements close to zero: ', num2str(nnz(ns(ns > eps) < 10. * TOL1))]);
disp(' ');

% Calculates restricted right null space to check behavior of clean
null_mat = null_mat(:, 1:7000);

tic;
sns = sp_null(null_mat, 1, TOL1);
T = toc;
err = max(max(abs(null_mat * sns)));
disp(['Found small null space in ', num2str(T), ...
      ' seconds with error of ', num2str(err)]);
disp(['Fractional fill is ', num2str(nnz(sns) / numel(sns))]);
disp(['Number of elements close to zero: ', num2str(nnz((sns < 10. * TOL1) & (sns > eps)))]);
disp(' ');

% Increasing size of null_mat can result in significant slowdowns due to
% memory requirements, even with fractional fill of around 0.025. This
% seems to be an unavoidable consequence of orthogonaliziation.
tic;
cns = clean(sns, TOL2);
T = toc;
err = max(max(abs(null_mat * cns)));
disp(['Found clean null space in ', num2str(T), ...
      ' seconds with error of ', num2str(err)]);
disp(['Fractional fill is ', num2str(nnz(cns) / numel(cns))]);
disp(['Number of elements close to zero: ', num2str(nnz((cns < 10. * TOL2) & (cns > eps)))]);