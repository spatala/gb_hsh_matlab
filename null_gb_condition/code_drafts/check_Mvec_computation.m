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

pts = 10;
c = zeros(pts*pts/2,1); 
p = zeros(pts*pts/2,1);
th1 = linspace(0,pi,pts/2); 
phi1 = linspace(0,2*pi,pts);
[phi1v, th1v] = meshgrid(phi1, th1);

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

ct1 = ceil(rand()*size(gb_rots,3));
gb = gb_rots(:,1:3,ct1);

s = convert_gbrots(gb_rots(:,:,ct1));
nsz = size(tot_inds,1);
Mvec0 = zeros(1,nsz);
for ct2=1:nsz
% ct2 = 5;

a  = tot_inds(ct2,3);
b  = tot_inds(ct2,4);
g1  = tot_inds(ct2,5);

al_val =  g1;
be_val = -g1;
alp_val = tot_inds(ct2,6);
bep_val =  tot_inds(ct2,7);


arr1 = tot_inds(ct2,3:7);


% a=1;
% b=0;
% g1=0;
% al_val =  g1;
% be_val = -g1;
% alp_val = -1;
% bep_val =  0;

a1_rots  = s.a1;  a2_rots  = s.a2;
lb1_rots = s.lb1; lb2_rots = s.lb2;
q1_rots  = s.q1;  q2_rots  = s.q2;
Q1_rots  = s.Q1;  Q2_rots  = s.Q2;

U1 = rotation_wo_svd(a1_rots,lb1_rots, q1_rots, Q1_rots, a, alp_val, al_val);
U2 = rotation_wo_svd(a2_rots,lb2_rots, q2_rots, Q2_rots, b, bep_val, be_val);

if a == b
    PI_ab = sqrt(2*a+1);
else
    if (a > b)
        PI_ab = sqrt(prod(2*b+1:2*a+1));
    else
        PI_ab = sqrt(prod(2*a+1:2*b+1));
    end
end

Mvec0(ct2) = PI_ab*(U1.*U2)/(sqrt(2*(pi^3)));
end



Mvec1 = calc_Mvec(gb,gb,symm_orders);

norm(Mvec0-Mvec1)