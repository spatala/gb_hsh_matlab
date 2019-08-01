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
symm_num = 14;
nrots = size(rots,3);

a_val = symm_orders(:,1)'; b_val = symm_orders(:,2)';
c_val = min(a_val, b_val);
num_cols = sum((2*a_val+1).*(2*b_val+1).*(2*c_val+1));


Mvecs = zeros(nrots,num_cols);

for ct1=1:nrots
    ct1
    g1_1 = rots(:,1:3,ct1); g1_2 = rots(:,4:6,ct1);
    ma1 = rots_to_angs(g1_1, g1_2);
    
    
    Mvec1 = zeros(1,num_cols);
    for a=a_val
        for b=b_val
            M1 = mbp_basis(a, b, [ma1(1), ma1(2), ma1(3), ma1(4), ma1(5)]);
            
            cond1 = tot_inds(:,3)==a & tot_inds(:,4)==b;
            ind_start = find(cond1,1);
            ind_stop  = find(cond1,1,'last');
            
            Mvec1(ind_start:ind_stop) = M1;
        end
    end
    Mvecs(ct1,:) = Mvec1;
end
Mvecs1 = sparse(Mvecs);

save(mat_name, 'Mvecs1', '-v7.3');
% symm_Mvec = Mvecs*S(:,symm_num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure('Position',[1200,50,1200,1200]); hold on;
% fval = abs(symm_Mvec);
% fval = reshape(fval, [num+1, num+1]);
% surf(X,Y,Z, fval);
% shading interp;
% axis equal; axis off;
% % view([1,0,0])
% view([1,1,1]);


rmpath(genpath(util_dir));