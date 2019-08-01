function mbp_angs = rots_to_angs(g1, g2)
mbp = rots_to_mbp([g1, g2]);
g_rots = mbp_to_rots(mbp);
g1 = g_rots(:,1:3); g2 = g_rots(:,4:6);
gb = g1; gm = (gb^(-1))*g2;

[w_m, th_m, ph_m] = rotmat_to_angs(gm);
[w_b, th_b, ph_b] = rotmat_to_angs(gb);
if (abs(th_b - pi/2) < 1e-14)
    mbp_angs = [w_m, th_m, ph_m, w_b, ph_b];
else
    error("theta_b is not equal to pi/2.")
end


end