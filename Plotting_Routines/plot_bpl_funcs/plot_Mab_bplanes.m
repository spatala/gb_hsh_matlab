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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num = 500;
[X, Y, Z] = sphere(num);
x1 = X(:); y1 = Y(:); z1 = Z(:);

% m1 = vrrotvec2mat([1,1,1, pi/3]);
% m1 = vrrotvec2mat([1,0,0, pi/4]);
m1 = vrrotvec2mat([1,0,0, 0]);
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

nrots = size(rots,3);
s = convert_gbrots(rots);

a1_rots  = s.a1;  a2_rots  = s.a2;
lb1_rots = s.lb1; lb2_rots = s.lb2;
q1_rots  = s.q1;  q2_rots  = s.q2;
Q1_rots  = s.Q1;  Q2_rots  = s.Q2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'Sarr_MabInds_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
tot_Uprops=s1.tot_Uprops; ind_ranges=s1.ind_ranges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_rows = size(S,1);


for ct1 = size(S,2)-1:size(S,2)
    ct1
% ct1 = 3;
st1 = ind_ranges(ct1,1);
st2 = ind_ranges(ct1,2);

U_props = tot_Uprops(st1:st2,:);
vec_inds = U_props(:,1);
a_val   = U_props(:,2); alp_val = U_props(:,5); al_val  = U_props(:,4);
b_val   =  U_props(:,3); bep_val =  U_props(:,6); be_val  = -U_props(:,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mvec = sparse(nrots,num_rows);
tic;
for ct2 = 1:size(vec_inds,1)
    U1 = rotation_wo_svd(a1_rots,lb1_rots, q1_rots, Q1_rots, a_val(ct2), alp_val(ct2), al_val(ct2));
    U2 = rotation_wo_svd(a2_rots,lb2_rots, q2_rots, Q2_rots, b_val(ct2), bep_val(ct2), be_val(ct2));
    Mvec(:,vec_inds(ct2)) = U1.*U2;
end
toc;


symm_Mvec = Mvec*S(:,ct1);

% ct1 = 9;
figure('Position',[1200,50,1200,1200]); hold on;

% X1 = pts(:,1)./(1+pts(:,3)); Y1 = pts(:,2)./(1+pts(:,3));
% tri = delaunay(X1,Y1);
% trisurf(tri, X1, Y1, abs(symm_Mvec(:,ct1)));
fval = abs(symm_Mvec);
fval = reshape(fval, [num+1, num+1]);
surf(X,Y,Z, fval);
shading interp;
axis equal; axis off;
view([1,0,0])
% view([1,1,1])

end
% % [th,ph] = meshgrid(th1,ph1);
% % 
% % x1 = sin(th).*cos(ph);
% % y1 = sin(th).*cos(ph);
% % z1 = cos(th);
% % 
% % surf(x1,y1,z1)

rmpath(genpath(util_dir));