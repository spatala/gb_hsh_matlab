function pts = gen_sph_pts(num)

% th1 = linspace(0,pi/2, num)';
% ph1 = linspace(0,2*pi, 4*num);

th1 = linspace(0,pi, num)';
ph1 = linspace(0,2*pi, 2*num);

x1 = sin(th1)*cos(ph1); x1 = x1(:);
y1 = sin(th1)*sin(ph1); y1 = y1(:);
z1 = cos(th1)*ones(size(ph1)); z1 = z1(:);


pts = [x1, y1, z1];
pts = unique(pts, 'rows');
%  

end