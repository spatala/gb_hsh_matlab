function [] = generate_gb_cryst_symm(pt_grp, Nmax, TOL)
% function [] = generate_gb_cryst_symm()
% pt_grp = 'Oh'; Nmax = 4; TOL = 1e-12;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%
disp('Step 1');
generate_cryst_symm_ab(top_dir, pt_grp, Nmax, TOL);
%%%%
disp('Step 2');
generate_ges_mat(top_dir, pt_grp, Nmax)
%%%%
disp('Step 3');
combine_cryst_ges(top_dir,pt_grp, Nmax, TOL)
%%%%
disp('Step 4');
generate_gb_null(top_dir,pt_grp, Nmax)
%%%%
disp('Step 5');
combine_cryst_ges_gbnull(top_dir,pt_grp, Nmax, TOL)
%%%%
disp('Step 6');
save_symmvec_MabInds(top_dir,pt_grp, Nmax)
%%%%
disp('Step 7');
rmpath(genpath(util_dir));
end

function [] = generate_ges_mat(top_dir, pt_grp, Nmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'aPLUSb_max_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_aPLUSb_max_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;

[~, ~, ~, Laue] = get_symmgen_mats(pt_grp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flip_mat = generate_flip_mat(symm_orders);
if Laue
    symm_mat = flip_mat;
else
    ypi_left = generate_ypi_left(symm_orders);
    symm_mat = ypi_left*flip_mat;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, 'ges_mat_aPLUSb_max_',num2str(Nmax),'.mat'];
save(mat_name,'symm_mat');
end

function combine_cryst_ges(top_dir,pt_grp, Nmax, TOL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'aPLUSb_max_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'ges_mat_aPLUSb_max_',num2str(Nmax),'.mat'];
s1 = load(mat_name); R1 = s1.symm_mat;

mat_name = [data_fname0, ...
    'Sarr_abc_combined_csymm_aPLUSb_max_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S0 = sparse(s1.S);

S = S0 * sp_orth(sp_null(S0' * R1 * S0 - eye(size(S0, 2)), 1, TOL));

mat_name = [data_fname0, ...
    'Sarr_cryst_ges_aPLUSb_max_',num2str(Nmax),'.mat'];
save(mat_name,'S');

end


function combine_cryst_ges_gbnull(top_dir,pt_grp, Nmax, TOL)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'aPLUSb_max_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'gbnull_mat_aPLUSb_max_',num2str(Nmax),'.mat'];
s1 = load(mat_name); R1 = s1.null_mat;

mat_name = [data_fname0, ...
    'Sarr_cryst_ges_aPLUSb_max_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S0 = sparse(s1.S);

S = S0 * sp_orth(sp_null(R1 * S0, 1, TOL));
tic; S = clean(S, TOL); toc;

disp(size(S));

mat_name = [data_fname0, ...
    'Sarr_cryst_ges_gbnull_aPLUSb_max_',num2str(Nmax),'.mat'];
save(mat_name,'S');

end



