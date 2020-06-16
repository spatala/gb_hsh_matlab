function rot_angs = rotmat_to_angs(g)
%
% Convert rotation maxtrix to angles
% 
% Input:
%    g: $3 \times 3$ rotation matrix.
% 
% Output:
%    rot_angs: $1 \times 3$ rotation angles
%        1) w: Omega (rotation angle)
%        2) th: polar angle ? (theta), and 
%        3) ph: azimuthal angle ?
%

axang = vrrotmat2vec(g);
[az, el, ~] = cart2sph(axang(1),axang(2),axang(3));
ph = az; th = pi/2 - el; w = axang(4);
rot_angs = [w, th, ph];
end