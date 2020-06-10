function [w, th, ph] = CK_to_angles(CK)
%
% returns the rotation angle and axis for the given U^{1/2} with rows 
% and cols ordered by increasing m' and m.
% 
% 
% Input
%    CK: Cayley-Klein Parameters
% 
% Output
%    w, th, ph:
%        The rotation angle, polar and azimuthal angles of the rotation
%        matrix.
%

    a = CK(2, 2);
    b = CK(2, 1);
    
    w  = 2. * acos(real(a));
    tmp = -sin(w / 2.);
    th = acos(imag(a) / tmp);
    tmp = tmp * sin(th);
    ph = atan2(real(b) / tmp, imag(b) / tmp);
end