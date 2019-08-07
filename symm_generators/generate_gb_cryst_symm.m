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
generate_cryst_symm_ab(top_dir, pt_grp, Nmax, TOL);
%%%%
generate_ges_mat(top_dir, pt_grp, Nmax)
%%%%
combine_cryst_ges(top_dir,pt_grp, Nmax, TOL)
%%%%
generate_gb_null(top_dir,pt_grp, Nmax)
%%%%
combine_cryst_ges_gbnull(top_dir,pt_grp, Nmax, TOL)
%%%%
save_symmvec_MabInds(top_dir,pt_grp, Nmax)
%%%%
rmpath(genpath(util_dir));
end



function [] = generate_ges_mat(top_dir, pt_grp, Nmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
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
mat_name = [data_fname0, 'ges_mat_nmax_',num2str(Nmax),'.mat'];
save(mat_name,'symm_mat');
end
 
function combine_cryst_ges(top_dir,pt_grp, Nmax, TOL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'ges_mat_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); R1 = s1.symm_mat;

mat_name = [data_fname0, ...
    'Sarr_abc_combined_csymm_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S0 = sparse(s1.S);

S = S0 * sp_orth(sp_null(S0' * R1 * S0 - eye(size(S0, 2)), 1, TOL));

mat_name = [data_fname0, ...
    'Sarr_cryst_ges_nmax_',num2str(Nmax),'.mat'];
save(mat_name,'S');

end


function combine_cryst_ges_gbnull(top_dir,pt_grp, Nmax, TOL)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'gbnull_mat_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); R1 = s1.null_mat;

mat_name = [data_fname0, ...
    'Sarr_cryst_ges_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S0 = sparse(s1.S);

S = S0 * sp_orth(sp_null(R1 * S0, 1, TOL));
tic; S = clean(S, TOL); toc;

disp(size(S));

mat_name = [data_fname0, ...
    'Sarr_cryst_ges_gbnull_nmax_',num2str(Nmax),'.mat'];
save(mat_name,'S');

end


function Ypi_symm_mat = generate_ypi_left(symm_orders)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1,b1);
num_cols = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
nsymm_ab = size(symm_orders,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ypi_symm_mat = sparse(num_cols, num_cols);

ind_start = 1;
for ct1 = 1:nsymm_ab
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a_val = symm_orders(ct1,1); b_val = symm_orders(ct1,2);
    Mrot = generate_ypi_left_ab(a_val,b_val);
    c_val = min(a_val, b_val);
    Na = 2*a_val; Nb = 2*b_val; Nc = 2*c_val; nsz = (Na+1)*(Nb+1);
    
    nsz1 = nsz*(Nc+1);
    ind_stop = ind_start + nsz1 - 1;
    Ypi_symm_mat(ind_start:ind_stop, ind_start:ind_stop) = Mrot;
    ind_start = ind_stop + 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end




function flip_mat = generate_flip_mat(symm_orders)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1 = symm_orders(:,1); b1 = symm_orders(:,2); c1 = min(a1,b1);
num_cols = sum((2*a1+1).*(2*b1+1).*(2*c1+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% The variable tot_inds contains the mapping between the column-index
%%%% for the M-function and the (a,b,gamma,alpha, beta) values.
tot_inds = mbp_inds_ab_array(symm_orders);

%%%% flip_mat is the $A$ matrix that specifies the grain-exchange-symmetry
%%%% for the coefficients.
flip_mat = sparse(num_cols, num_cols);
% ind_maps = zeros(num_cols,2);
for ct1=1:num_cols
    %%%% Get the (a,b,gamma,alpha, beta) corresponding to the column.
    a1 = tot_inds(ct1,3); b1 = tot_inds(ct1,4);
    gamma1 = tot_inds(ct1,5); 
    alpha1 = tot_inds(ct1,6); beta1 = tot_inds(ct1,7);
    %%%% Find the index for (b,a,-gamma,beta, alpha) row.
    a2 = b1; b2 = a1; gamma2 = -gamma1;
    alpha2 = beta1; beta2 = alpha1;
    ind1 = find(...
        (tot_inds(:,3) == a2     ) & ...
        (tot_inds(:,4) == b2     ) & ...
        (tot_inds(:,5) == gamma2 ) & ...
        (tot_inds(:,6) == alpha2 ) & ...
        (tot_inds(:,7) == beta2  ));
    
    %%%% The (a,b, gamma, alpha, beta)th column contains a 1 in the 
    %%%% (b,a, -gamma, beta, alpha)th row
    flip_mat(ind1,ct1) = 1;
%     ind_maps(ct1,:) = [ct1, ind1];
end

end



