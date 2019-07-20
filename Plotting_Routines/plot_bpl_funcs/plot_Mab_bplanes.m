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

pt_grp = 'Oh'; Nmax = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1 = vrrotvec2mat([1,1,1, pi/3]);
% m1 = vrrotvec2mat([1,0,0, pi/4]);
% m1 = vrrotvec2mat([1,0,0, 0]);
% q1 = [cos(pi/4), 0, 0, sin(pi/4)];
% m1 = quat2mat(q1);
num = 500;
[X, Y, Z, rots] = gen_gb_bpl_rots(num, m1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
mat_name = [data_fname0, 'Sarr_cryst_ges_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
%%%% 
symm_num = 9;
symm_Mvec = compute_symm_Mvec(rots, S, symm_num, data_fname0, Nmax);
%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',[1200,50,1200,1200]); hold on;
fval = abs(symm_Mvec);
fval = reshape(fval, [num+1, num+1]);
surf(X,Y,Z, fval);
shading interp;
axis equal; axis off;
% view([1,0,0])
view([1,1,1])


rmpath(genpath(util_dir));