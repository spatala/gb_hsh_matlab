function mbp_angs = rots_to_angs(g1, g2)
mbp = rots_to_mbp([g1, g2]);
g_rots = mbp_to_rots(mbp);
g1 = g_rots(:,1:3); g2 = g_rots(:,4:6);
gb = g1; gm = (gb^(-1))*g2;

ra_m = rotmat_to_angs(gm);
ra_b = rotmat_to_angs(gb);
if (abs(ra_b(2) - pi/2) < 1e-14)
    mbp_angs = [ra_m, ra_b(1), ra_b(3)];
else
    error("theta_b is not equal to pi/2.")
end


end