function [Ylm] = spherical_harmonics(l, th, ph)
    Ylm = zeros(2 * l + 1, 1);
    
    th2 = mod(th, pi);
    ph2 = mod(ph, 2. * pi);
    
    % Construct associated Legendre functions of all orders
    tmp = legendre(l, cos(th2));
    Plm = [flipud(tmp(2:(l + 1))); tmp];
    for m = -l:-1
        tmp = ((-1)^(-m)) * (prod(1:(l + m)) / prod(1:(l - m)));
        Plm(m + l + 1) = tmp * Plm(m + l + 1);
    end
    
    % Construct spherical harmonics
    for m = -l:l
        tmp = sqrt(prod(1:(l - m)) / prod(1:(l + m))); 
        Ylm(l + m + 1) = tmp * Plm(m + l + 1) * exp(1i * m * ph2);
    end
    
    tmp = sqrt((2 * l + 1) / (4 * pi));
    Ylm = tmp * Ylm;
    
    % Phase change when th outside allowed bounds
    n = round(abs(th - th2) / pi);
    if mod(n, 2) == 1
        Ylm = (-1)^l * Ylm;
    end
    
    % Phase change when ph outside allowed bounds
    n = round(abs(ph - ph2) / pi);
    if mod(n, 2) == 1
        Ylm = (-1)^(-l:l) .* Ylm;
    end
end