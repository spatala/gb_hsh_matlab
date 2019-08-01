function check_conv(r1, rot_angs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Check angles match with rotation matrix
w = rot_angs(1); th = rot_angs(2); ph = rot_angs(3);

ax1 = [sin(th)*cos(ph),sin(th)*sin(ph),cos(th)];
if (norm(r1 - vrrotvec2mat([ax1,w])) > 1e-12)
    error("Conversion from rotation matrix to polar angles failed.")
end

end