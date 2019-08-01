function [w, th, ph] = rotmat_to_angs(g)
axang = vrrotmat2vec(g);
[az, el, ~] = cart2sph(axang(1),axang(2),axang(3));
ph = az; th = pi/2 - el; w = axang(4);
end