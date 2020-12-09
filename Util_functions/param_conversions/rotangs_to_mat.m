function g = rotangs_to_mat(angs)
% 
% Convert angles to rotation maxtrix
% 
% - Input:
%   + rot_angs: $1 \times 3$ rotation angles
% 	- w: Omega (rotation angle)
% 	- th: polar angle (theta), and 
% 	- ph: azimuthal angle (phi)
% 
% - Output:
%   + g: $3 \times 3$ rotation matrix.
% 

w = angs(1); th = angs(2); ph = angs(3);

z1 = cos(th);
x1 = sin(th)*cos(ph);
y1 = sin(th)*sin(ph);

g = vrrotvec2mat([x1,y1,z1,w]);

end