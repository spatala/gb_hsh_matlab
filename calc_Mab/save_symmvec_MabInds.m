clear all; clc;

pt_grp = 'Oh'; Nmax = 10;

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

% s1 = set_vars(); Nmax = s1.Nmax; pt_grp = s1.pt_grp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'nmax_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_Nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_nmax_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
S = s1.S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Move this to a separate code.
num_rows = size(S,1);
num_cols = size(S,2);
tot_inds = mbp_inds_ab_array(symm_orders, num_rows);
tot_Uprops = zeros(num_rows*num_cols,6);
ind_ranges = zeros(num_cols,2);
st1 = 1;
for ct1 = 1:num_cols
    ct1
    ind1 = find(abs(S(:,ct1))>0);
    st2 = st1 - 1 + size(ind1,1);
    U_props = tot_inds(ind1,:);
    tot_Uprops(st1:st2,:) = U_props(:,[1,3:end]);
    ind_ranges(ct1,:) = [st1, st2];
    st1 = st2 + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
tot_Uprops(st2+1:end,:) = [];
[C, ia, ic] = unique(tot_Uprops, 'rows');
%%% C = A(ia,:) and A = C(ic,:)

mat_name = [data_fname0, ...
    'Sarr_MabInds_nmax_',num2str(Nmax),'.mat'];
save(mat_name, 'tot_Uprops', 'ind_ranges');


% a_val   = U_props(:,3); alp_val = U_props(:,6); al_val  = U_props(:,5);
% b_val   =  U_props(:,4); bep_val =  U_props(:,7); be_val  = -U_props(:,5);
% rmpath(genpath(util_dir));