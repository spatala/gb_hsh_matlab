function mbp = rots_to_mbp(rots)
%
% Converts GB parameters from SO3xSO3 parameterization (rots) to 
% Misorientation-Boundary plane (MBP).
% 
% - Input
% 	+ rots: 1 X 3 X 6 matrix
%       - First three columns are the first orientation matrix (g1) and
%       - the second three are the second orientaiton matrix (g2)
% 
% - Output
% 	+ mbp: 3 X 4 matrix
%       - First three columns are the misorientation matrix and 
%       - the last column is the boundary-plane.
%

g1 = rots(:,1:3); g2 = rots(:,4:6);

g1_inv = g1^(-1);
zvec = [0,0,1]';
m1 = g1_inv*g2;
n1 = g1_inv*zvec;
mbp = [m1, n1];
end