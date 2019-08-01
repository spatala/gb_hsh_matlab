function [CK] = angles_to_CK(w, th, ph)
% [CK] = angles_to_CK(w, th, ph) - returns U^{1/2} for the rotation given
%   by a rotation angle and axis. Rows and cols ordered by increasing m'
%   and m.
    a = cos(w / 2.) - 1i * sin(w / 2.) * cos(th);
    b = -1i * sin(w / 2.) * sin(th) * exp(-1i * ph);
    
    CK = [conj(a), -conj(b); b, a];
end