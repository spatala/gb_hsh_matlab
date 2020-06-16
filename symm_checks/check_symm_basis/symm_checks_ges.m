function [] = symm_checks_ges(pt_grp, Nmax, coeffs_typ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Computes the full MBP basis function for a random GB parameter 
% 2) Computes the full MBP basis function for a GB parameter with GES
% 3) Checks that the GES symmetrized basis function for the two GB
% parameters are equivalent.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
top_dir = get_top_dir();
util_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(util_dir));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Get data files
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,coeffs_typ,'_',num2str(Nmax),'/'];

%%%%% Load the possible (a,b) values for (pt_grp, Nmax)
mat_name = [data_fname0,'symm_ab_',pt_grp,'_',coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
a_val = symm_orders(:,1)'; b_val = symm_orders(:,2)'; c_val = min(a_val, b_val);
num_cols = sum((2*a_val+1).*(2*b_val+1).*(2*c_val+1));
%%%% tot_inds gives a mapping between (a,b) and range of indices in the
%%%% list of basis functions.
%%%% This should be made faster!
tot_inds = mbp_inds_ab_array(symm_orders);

%%%% Load mat file containing the symmetrized basis functions
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_gbnull_',coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
nsymm_evs = size(S,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Check GES "num" GB parameters
num = 100;
diff_norms = zeros(num,1);
for tct1 = 1:num
    %%%% Generate Random Grain boundary (g1_1, g1_2)
    diff_vec = zeros(nsymm_evs,1);
    w1 = rand()*pi; th1 = rand()*pi; ph1 = 2. * pi * rand();
    axang1 = [sin(th1)*cos(ph1), sin(th1)*sin(ph1), cos(th1), w1];
    g1_1 = vrrotvec2mat(axang1);
    w2 = rand()*pi; th2 = rand()*pi; ph2 = 2. * pi * rand();
    axang2 = [sin(th2)*cos(ph2), sin(th2)*sin(ph2), cos(th2), w2];
    g1_2 = vrrotvec2mat(axang2);
    ma1 = rots_to_angs(g1_1, g1_2);

    %%%% Generate GES of Random Grain boundary (g2_1, g2_2)
    Ypi = vrrotvec2mat([0,1,0,pi]);
    g2_1 = Ypi*g1_2; g2_2 = Ypi*g1_1;
    ma2 = rots_to_angs(g2_1, g2_2);
    
    %%%% Compute Mvec1 and Mvec2
    Mvec1 = zeros(1,num_cols); Mvec2 = zeros(1,num_cols);
    for ct2 = 1:size(symm_orders,1)
        a = a_val(ct2); b = b_val(ct2);
        M1 = mbp_basis(a, b, ma1);
        M2 = mbp_basis(a, b, ma2);
        
        cond1 = tot_inds(:,3)==a & tot_inds(:,4)==b;
        ind_start = find(cond1,1);
        ind_stop  = find(cond1,1,'last');
        
        Mvec1(ind_start:ind_stop) = M1;
        Mvec2(ind_start:ind_stop) = M2;
    end
    
    %%%%% Compute the difference between two Mvectors
    for ct1=1:nsymm_evs
        diff_vec(ct1) = norm(Mvec1*S(:,ct1) - Mvec2*S(:,ct1));
    end
    diff_norms(tct1) = norm(diff_vec)/nsymm_evs;
end

norm(diff_norms)/num

rmpath(genpath(util_dir));