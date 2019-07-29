clear all; clc;

th = rand()*pi; phi = rand()*2*pi;
l = ceil(rand()*10);


Yl = compute_complex_Yl(th,phi,l);

% Yl_vals = zeros(2*l+1,1);
% ct1 = 1;
% for m = -l:l
%     Y1 = harmonicY(l,m,th,phi,'type','complex', 'norm', true, 'phase', true);
%     Yl_vals(ct1) = Y1;
%     ct1 = ct1 + 1;
% end
% 
% norm(Yl-Yl_vals)

%%%% Compute multiplication factor
a1 = 3; b1 = 2; al1 = 2; be1 = -1; g1 = 1;

