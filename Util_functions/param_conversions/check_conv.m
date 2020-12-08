function check_conv(r1, rot_angs)
% 
% Function checks that the angles (w, th, ph) match with the rotation
% matrix (r1)
% 
% - Input:
%   + r1: Rotation matrix
%   + rot_angs: array of rotation angles (w, th, ph)
% 
% - Output: None
%   + Prints an error message if the check fails.
% 
w = rot_angs(1); th = rot_angs(2); ph = rot_angs(3);

ax1 = [sin(th)*cos(ph),sin(th)*sin(ph),cos(th)];
if (norm(r1 - vrrotvec2mat([ax1,w])) > 1e-10)
    error("Conversion from rotation matrix to polar angles failed.")
end

end