function Yl_vals = compute_complex_Yl(th,phi,l)

% include Condon-Shortley phase term
% this term is in fact included in MATLAB's associated Legendre function
% definition, so what we are doing here is removing it if not requested
CS_phase = true;

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
    
    % normalization
    b = 1/prod(l-m+1:l+m); % factorial(n-m)/factorial(n+m);
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
end
