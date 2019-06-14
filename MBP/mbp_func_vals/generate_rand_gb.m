clear all; clc;

num_pts = 2;


sph_pts = zeros(num_pts,3);
u1 = 2*rand(num_pts, 1)-1;
th = 2*pi*rand(num_pts,1);
sph_pts(:,1) = sqrt(1-u1.^2).*cos(th);
sph_pts(:,2) = sqrt(1-u1.^2).*sin(th);
sph_pts(:,3) = u1;


quat_pts = zeros(num_pts,4);
u_pts = rand(num_pts, 3);
u1 = u_pts(:,1);u2 = u_pts(:,2);u3 = u_pts(:,3);
quat_pts(:,1) = sqrt(1-u1).*sin(2*pi*u2);
quat_pts(:,2) = sqrt(1-u1).*cos(2*pi*u2);
quat_pts(:,3) = sqrt(u1).*sin(2*pi*u3);
quat_pts(:,4) = sqrt(u1).*cos(2*pi*u3);

gbs = [quat_pts, sph_pts];