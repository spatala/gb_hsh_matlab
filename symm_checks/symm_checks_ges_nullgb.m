%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Code to check the rotation of SO(4) Functions R^(a,b)
%%%%%%%
clear all; clc;

% pt_grp = 'C1'; Nmax = 1;
pt_grp = 'Oh'; Nmax = 4;

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

% s1 = set_vars();
% Nmax = s1.Nmax; pt_grp = s1.pt_grp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tot_inds = mbp_inds_ab_array(symm_orders);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_gbnull_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
nsymm_evs = size(S,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Check Null-GB condition for zero misorientation
w1 = rand()*pi; th1 = rand()*pi; ph1 = 2. * pi * rand();
axang1 = [sin(th1)*cos(ph1), sin(th1)*sin(ph1), cos(th1), w1];
g1_1 = vrrotvec2mat(axang1); g1_2 = g1_1;
ma1 = rots_to_angs(g1_1, g1_2);
% mbp_angs(tct1,1:5) = ma1;
    
a_val = symm_orders(:,1)';
b_val = symm_orders(:,2)';
c_val = min(a_val, b_val);
num_cols = sum((2*a_val+1).*(2*b_val+1).*(2*c_val+1));

diff_vec = zeros(nsymm_evs,1);
for ct1=1:nsymm_evs

Mvec = zeros(1,num_cols);

for a=a_val
    for b=b_val
        M1 = mbp_basis(a, b, ma1(1), ma1(2), ma1(3), ma1(4), ma1(5));

        cond1 = tot_inds(:,3)==a & tot_inds(:,4)==b;
        ind_start = find(cond1,1);
        ind_stop  = find(cond1,1,'last');

        Mvec(ind_start:ind_stop) = M1;
    end
end
    diff_vec(ct1) = norm(Mvec*S(:,ct1));
end
norm(diff_vec)

num = 10;
diff_norms = zeros(num,1);
mbp_angs = zeros(num,10);
for tct1 = 1:num
    tct1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Check GES for finite misorientation
    diff_vec = zeros(nsymm_evs,1);
    w1 = rand()*pi; th1 = rand()*pi; ph1 = 2. * pi * rand();
    axang1 = [sin(th1)*cos(ph1), sin(th1)*sin(ph1), cos(th1), w1];
    g1_1 = vrrotvec2mat(axang1);
    w2 = rand()*pi; th2 = rand()*pi; ph2 = 2. * pi * rand();
    axang2 = [sin(th2)*cos(ph2), sin(th2)*sin(ph2), cos(th2), w2];
    g1_2 = vrrotvec2mat(axang1);
    ma1 = rots_to_angs(g1_1, g1_2);
    mbp_angs(tct1,1:5) = ma1;
    
    
    Ypi = vrrotvec2mat([0,1,0,pi]);
    g2_1 = Ypi*g1_2; g2_2 = Ypi*g1_1;
    ma2 = rots_to_angs(g2_1, g2_2);
    mbp_angs(tct1,6:10) = ma2;
    
    a_val = symm_orders(:,1)';
    b_val = symm_orders(:,2)';
    c_val = min(a_val, b_val);
    num_cols = sum((2*a_val+1).*(2*b_val+1).*(2*c_val+1));
    
    mbp1 = rots_to_mbp([g1_1,g1_2]);
    mbp2 = rots_to_mbp([g2_1,g2_2]);
    
    M1 = mbp1(:,1:3); n1 = mbp1(:,4);
    M2 = mbp2(:,1:3); n2 = mbp2(:,4);
    
    % [norm(M1 - M2'), norm(n2 + M1'*n1)]
    
    for ct1=1:nsymm_evs
        Mvec1 = zeros(1,num_cols); Mvec2 = zeros(1,num_cols);
        for a=a_val
            for b=b_val
                M1 = mbp_basis(a, b, ma1(1), ma1(2), ma1(3), ma1(4), ma1(5));
                M2 = mbp_basis(a, b, ma2(1), ma2(2), ma2(3), ma2(4), ma2(5));
                
                cond1 = tot_inds(:,3)==a & tot_inds(:,4)==b;
                ind_start = find(cond1,1);
                ind_stop  = find(cond1,1,'last');
                
                Mvec1(ind_start:ind_stop) = M1;
                Mvec2(ind_start:ind_stop) = M2;
            end
        end
        diff_vec(ct1) = norm(Mvec1*S(:,ct1) - Mvec2*S(:,ct1));
    end
    diff_norms(tct1) = norm(diff_vec);
end

rmpath(genpath(util_dir));