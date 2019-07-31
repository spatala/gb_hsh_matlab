clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% addpath(genpath('Spherical-Harmonic-Transform/'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pt_grp = 'C1'; Nmax = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tot_inds = mbp_inds_ab_array(symm_orders);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Get boundary-planes for Identity misorientation
pts = 10;
c = zeros(pts*pts/2,1); 
p = zeros(pts*pts/2,1);
th1 = linspace(0,pi,pts/2); 
phi1 = linspace(0,2*pi,pts);
[phi1v, th1v] = meshgrid(phi1, th1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Conver (I, n) to (g1, g1)
%%%% Boundary normals are given by phi2, th2
phi2 = phi1v(:); th2 = th1v(:);
bpn_vecs = [sin(th2).*cos(phi2), sin(th2).*sin(phi2), cos(th2)];
[bpn_vecs, ia, ic] = uniquetol(bpn_vecs,1e-10,'ByRows',true);
phi2 = phi2(ia); th2 = th2(ia);

nvecs = size(bpn_vecs,1);
gb_rots = zeros(3,6,nvecs);
om_b = zeros(nvecs,1); phi_b = zeros(nvecs,1);
for ct1 = 1:nvecs
    mbp = [eye(3), bpn_vecs(ct1,:)'];
    gb_rots(:,:,ct1) = mbp_to_rots(mbp);
    ax_ang = vrrotmat2vec(gb_rots(:,1:3,ct1));
    ax1 = ax_ang(1:3);
    om_b(ct1) = ax_ang(4);
    [phi_b(ct1), ~, ~] = cart2pol(ax1(1), ax1(2), ax1(3));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ct1 = 2; %%% Index for boundary-plane iteration
gb = gb_rots(:,1:3,ct1); om_b1 = om_b(ct1); phi_b1 = phi_b(ct1);
%%%% Compute Mvec using SO(3) irreps
Mvec = calc_Mvec(gb,gb,symm_orders);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ct2 = 2; %%% Index for (a,b) iteration
a1  = tot_inds(ct2,3); b1  = tot_inds(ct2,4); arr1 = tot_inds(ct2,3:7);
Z = 0;
for e_val = abs(a1-b1):(a1+b1)
    Z = Z + (sph_harm_Mab_formula(arr1,e_val,om_b1, phi_b1));
end
% PI_ab = get_PIab(a1,b1);
PI_ab = (sqrt(2*a1+1))*(sqrt(2*b1+1));
Z = sqrt(2)*PI_ab*Z/pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

norm(Z-Mvec(ct2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rmpath(genpath(util_dir));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ct_f = 0;
% for ct2=1:size(tot_inds,1)
% a = 3; b = 2; g1 = 2; al1 = 3; be1 = -1;
% g1  = tot_inds(ct2,5);
% al1 = tot_inds(ct2,6);
% be1 = tot_inds(ct2,7);




% if (abs(Z-Mvec(ct2)) < 1e-13)
%     display(arr1);
% else
%     ct_f = ct_f + 1;
% end
% 
% end

