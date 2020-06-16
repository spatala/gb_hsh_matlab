function mbp_angs = rots_to_angs(g1, g2)
%
% Convert $SO(3) \times SO(3)$ parameterization into five angles
% 
% - Input
%   + g1, g2: $3 \times 3$ rotation matrices.
% 
% - Output
%   + mbp_angs: 1 X 5 matrix. The columns are as follows:
% 	- omega_m: rotation angle
% 	- theta_m: polar angle
% 	- phi_m  : azimuthal angle
% 	- omega_b: rotation angle
% 	- phi_b  : azimuthal angle
% 
% - Notes
% 	+ The first three angles give the misorientation matrix. 
%   The last two give the rotation matrix to bring the boundary-plane 
%   normal parallel to z-axis. The rotation is along a vector that 
%   is in the X-Y plane.
%

% g1 = rots(:,1:3); g2 = rots(:,4:6);
% mbp = rots_to_mbp([g1, g2]);
% g_rots = mbp_to_rots(mbp);
% g1 = g_rots(:,1:3); g2 = g_rots(:,4:6);
gb = g1; gm = (gb^(-1))*g2;

ra_m = rotmat_to_angs(gm);
ra_b = rotmat_to_angs(gb);
if (abs(ra_b(2) - pi/2) < 1e-14)
    mbp_angs = [ra_m, ra_b(1), ra_b(3)];
else
    error("theta_b is not equal to pi/2.")
end


end