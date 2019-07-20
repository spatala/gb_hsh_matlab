function [X, Y, Z, rots] = gen_gb_bpl_rots(num, m1)

[X, Y, Z] = sphere(num);
x1 = X(:); y1 = Y(:); z1 = Z(:);

nsz = (num+1)^2;
rots = zeros(3,6,nsz);
for ct1 = 1:nsz
    n1 = [x1(ct1); y1(ct1); z1(ct1)];
    tmbp = [m1, n1];
    trots = mbp_to_rots(tmbp);
    rots(:,:,ct1) = trots;
end

% % th1 = linspace(0,pi/2, num)';
% % ph1 = linspace(0,2*pi, 4*num);
% 
% th1 = linspace(0,pi, num)';
% ph1 = linspace(0,2*pi, 2*num);
% 
% x1 = sin(th1)*cos(ph1); x1 = x1(:);
% y1 = sin(th1)*sin(ph1); y1 = y1(:);
% z1 = cos(th1)*ones(size(ph1)); z1 = z1(:);
% 
% 
% pts = [x1, y1, z1];
% pts = unique(pts, 'rows');
% %  

end