function [] = symm_checks_nullgb_constraint(pt_grp, Nmax, coeffs_typ, constraint)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Computes the full MBP basis function for a random GB parameter with
% identity misorientation
% 2) Computes the norm of the "constrained" basis-function
% 3) Checks that the norm is equal to zero.
%       Null boundary singularity where f(g,g) = 0;
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
tot_inds = mbp_inds_ab_array(symm_orders);

%%%% Load mat file containing the symmetrized basis functions
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_gbnull_',constraint,'_',...
    coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
nsymm_evs = size(S,2);

num1 = 1000;
null_rots = zeros(3,6,num1);

for ct1 = 1:num1
    %%%% Set-up Identity misorientation
    w_m  = 0; th_m = pi * rand(); ph_m = 2. * pi * rand();
    axang_m = [sin(th_m)*cos(ph_m), sin(th_m)*sin(ph_m), cos(th_m), w_m];
    gm = vrrotvec2mat(axang_m);
    
    %%%% Set-up random boundary-plane orientations
    w_b  = 8 * pi * (rand() - 1); th_b = pi / 2; ph_b = 2. * pi * rand();
    axang_b = [sin(th_b)*cos(ph_b), sin(th_b)*sin(ph_b), cos(th_b), w_b];
    gb = vrrotvec2mat(axang_b);
    g1 = gb;
    g2 = gb*gm;
    
    %     g1 = rotangs_to_mat([w_m, th_m, ph_m]);
    %     g2 = rotangs_to_mat([w_b, th_b, ph_b]);
    null_rots(:,:,ct1) = [g1, g2];
end

s_angs = convert_gbrots(null_rots);

% constraint = 'const';
SMvec=calc_Mvec_symm(top_dir, pt_grp, Nmax, ...
    coeffs_typ, constraint, s_angs);

if strcmp(constraint,'const')
    nsz = num1*(num1+1)/2;
    diff_vec = zeros(nsz,1);
    tct1 = 1;
    for ct1 = 1:num1-1
        % disp(ct1);
        for ct2 = ct1+1:num1
            diff_vec(tct1) = norm(full(SMvec(ct1,:)-SMvec(ct2,:)));
            tct1 = tct1 + 1;
        end
    end
    disp(max(diff_vec));
end

if strcmp(constraint,'zero')
    diff_vec = zeros(num1,1);
    for ct1 = 1:num1-1
        diff_vec(ct1) = norm(full(SMvec(ct1,:)));
    end
    disp(max(diff_vec));
end
rmpath(genpath(util_dir));
end