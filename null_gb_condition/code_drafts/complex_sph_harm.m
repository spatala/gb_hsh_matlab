function Z = complex_sph_harm(th2, phi2)
%%%% 
% clear all; clc;
% L_max = 8;
% pts = 10;
% c = zeros(pts*pts/2,1); 
% p = zeros(pts*pts/2,1);
% phi1 = linspace(0,pi,pts/2); 
% th1 = linspace(0,2*pi,pts);
% [phi1v, th1v] = meshgrid(phi1, th1);
% phi2 = phi1v(:); th2 = th1v(:);


for L = L_max
    Z = zeros(pts*pts/2, 1);
    P = legendre(L,cos(th2));
    for M = -L:L
        p = P(abs(M)+1,:)';
        Z = Z + realsqrt((2*L+1)/(4*pi*prod((L-M+1):(L+M))))*(p.*(cos(M*phi2) + 1i*sin(M*phi2)));
    end
end
end