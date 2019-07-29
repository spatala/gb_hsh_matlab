function Yl_vals = compute_complex_Yl(th,phi,l)
% clear all; clc;

% th = rand()*pi; phi = rand()*2*pi;
% th = 1.160022974928946;
% phi = 0.698707518179812;

CS_phase = true;

% l = ceil(rand()*10);
% m_vals = -l:l;
% m = m_vals(ceil(rand()*(2*l+1)));
% l = 3; m = 1;
% Y1 = harmonicY(l,m,th,phi,'type','complex', 'norm', true, 'phase', true);

%%%%% Dependent on $l$
a = (2*l+1)/(4*pi);
% associated Legendre function
P = legendre(l,cos(th));

Yl_vals = zeros(2*l+1,1);

ct1 = 1;
for m = -l:l
%%%%% Dependent on $m$
% check if m is odd and negative
isoddm = mod(m,2) == 1;
isnegm = m < 0;
% if m is negative, compute the symmetric function where m > 0
if isnegm
    m = abs(m);
end

% b = factorial(n-m)/factorial(n+m);
b = 1/prod(l-m+1:l+m);
C = sqrt(a*b);

Pm = P(abs(m)+1,:)';
E = exp(1i*m*phi);

% surface spherical harmonics
Y = C * Pm .* E;

% include Condon-Shortley phase term
% this term is in fact included in MATLAB's associated Legendre function
% definition, so what we are doing here is removing it if not requested
if ~CS_phase && isoddm
    Y = -Y;
end

% if m was negative
if isnegm
    Y = conj(Y);
    if isoddm
        Y = -Y;
    end
end
Yl_vals(ct1) = Y;
ct1 = ct1 + 1;
end

% [l, m, E, Y]
% Y - Y1

end
