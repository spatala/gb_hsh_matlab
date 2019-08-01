function rot_angs = rotmat_to_angs(g)
axang = vrrotmat2vec(g);
[az, el, ~] = cart2sph(axang(1),axang(2),axang(3));
ph = az; th = pi/2 - el; w = axang(4);
rot_angs = [w, th, ph];
end