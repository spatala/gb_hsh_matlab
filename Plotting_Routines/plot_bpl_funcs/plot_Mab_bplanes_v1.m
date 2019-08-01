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

pt_grp = 'Oh'; Nmax = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1 = vrrotvec2mat([1,1,1, pi/3]);
% m1 = vrrotvec2mat([1,0,0, pi/4]);
% m1 = vrrotvec2mat([1,0,0, 0]);
% q1 = [cos(pi/4), 0, 0, sin(pi/4)];
% m1 = quat2mat(q1);
num = 250;
[X, Y, Z, rots] = gen_gb_bpl_rots(num, m1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
mat_name = [data_fname0, 'Sarr_cryst_ges_gbnull_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
tot_inds = mbp_inds_ab_array(symm_orders);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mat_name = ['Mvecs_',num2str(num),'.mat'];


s1 = load(mat_name);
Mvecs=s1.Mvecs1;

symm_num = 1;
symm_Mvec = Mvecs*S(:,symm_num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',[1200,50,1200,1200]); hold on;
fval = abs(symm_Mvec);
fval = reshape(fval, [num+1, num+1]);
surf(X,Y,Z, fval);
shading interp;
axis equal; axis off;
% view([1,0,0])
view([1,1,1]);
colorbar;

rmpath(genpath(util_dir));