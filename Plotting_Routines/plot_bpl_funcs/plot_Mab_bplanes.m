clear all; clc;

curr_pwd = split(pwd,'/');
top_dir = '';
for ct1=1:length(curr_pwd)
    top_dir = strcat(top_dir,curr_pwd{ct1},'/');
    if (strcmp(curr_pwd{ct1},'gb_hsh_matlab'))
        break;
    end
end
util_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(util_dir));

s1 = set_vars(); Nmax = s1.Nmax; pt_grp = s1.pt_grp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% pts = gen_sph_pts(25);

% %%%%%% Rotate points
% %%% Z-axis goes to [1,1,1]
% zrot = [1,1,1]'; zrot = zrot/norm(zrot);
% %%% X-axis goes to [2,-1,-1]
% xrot = [2,-1,-1]'; xrot = xrot/norm(xrot);
% yrot = cross(zrot, xrot);
% 
% grot = [xrot, yrot, zrot];
% 
% rot_pts = (grot*pts')';
% % x1 = rot_pts(:,1); y1 = rot_pts(:,1); z1 = rot_pts(:,1);
% x1 = rot_pts(:,1); y1 = rot_pts(:,1); z1 = rot_pts(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num = 100;
[X, Y, Z] = sphere(num);
x1 = X(:); y1 = Y(:); z1 = Z(:);

m1 = vrrotvec2mat([1,1,1, pi/3]);
% m1 = vrrotvec2mat([1,0,0, pi/4]);
% m1 = vrrotvec2mat([1,0,0, 0]);
% q1 = [cos(pi/4), 0, 0, sin(pi/4)];
% m1 = quat2mat(q1);
nsz = (num+1)^2;
rots = zeros(3,6,nsz);
for ct1 = 1:nsz
    n1 = [x1(ct1); y1(ct1); z1(ct1)];
    tmbp = [m1, n1];
    trots = mbp_to_rots(tmbp);
%     mbp = rots_to_mbp(trots);
%     norm(mbp - tmbp)
    rots(:,:,ct1) = trots;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
num_cols = size(S,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mvec = zeros(nsz, num_cols);
parfor ct1=1:nsz
%     ct1
    rots1 = rots(:,:,ct1);
    g1 = rots1(:,1:3); g2 = rots1(:,4:6);
    tMvec = calc_Mvec(g1,g2,symm_orders);
    Mvec(ct1,:) = tMvec;
end

symm_Mvec = Mvec*S;

ct1 = 6;
figure(); hold on;

% X1 = pts(:,1)./(1+pts(:,3)); Y1 = pts(:,2)./(1+pts(:,3));
% tri = delaunay(X1,Y1);
% trisurf(tri, X1, Y1, abs(symm_Mvec(:,ct1)));
fval = abs(symm_Mvec(:,ct1));
fval = reshape(fval, [num+1, num+1]);
surf(X,Y,Z, fval);
shading interp;
axis equal; axis off;
view([1,0,0])

% % [th,ph] = meshgrid(th1,ph1);
% % 
% % x1 = sin(th).*cos(ph);
% % y1 = sin(th).*cos(ph);
% % z1 = cos(th);
% % 
% % surf(x1,y1,z1)

rmpath(genpath(util_dir));