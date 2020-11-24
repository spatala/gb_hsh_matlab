function save_symmvec_MabInds_gbnull_zero(top_dir, pt_grp, Nmax)
% 
% Save, "tot_Uprops" and "ind_ranges", that store the indices corresponding
% to non-zero components of the symmetrized basis functions.
% 
% - Input:
%   + top_dir: Path to identify data_files that contains symmetrized basis
%   functions.
%   + pt_grp: Point group symmetry of the underlying crystal.
%   + Nmax: Order of the basis functions (a,b) s.t. max(a+b) \leq Nmax.
% 
% - Output (stored in .mat file):
%   + tot_Uprops: Indices corresponding to the non-zero components of the
%   basis-vector. This array contains all the non-zero components of all 
%   the symmetrized basis-vector. They are vertically stacked. The 
%   index-range for each basis-vector is provided in ind_ranges variable.
%   The columns are 1(#),3(a),4(b),5(gamma),6(alpha),7(beta) of tot_inds.
% 
%   + ind_ranges: The number of rows are the number of symmetrized basis
%   functions. Each row corresponds to the "start" and "stop" indices in 
%   the tot_Uprops array containing the non-zero components.
%   For example, if the third-row in "ind_ranges" is [296, 544]. Then the
%   rows 296:544 contain the (a,b,gamma,alpha,beta) indices of non-zero
%   components for the "third" basis vector.
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,'aPLUSb_max_',num2str(Nmax),'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0,'symm_ab_',pt_grp,'_aPLUSb_max_',...
    num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mat_name = [data_fname0, ...
    'Sarr_cryst_ges_gbnull_zero_aPLUSb_max_',num2str(Nmax),'.mat'];
s1 = load(mat_name);
S = s1.S;

% num_rows: Number of basis functions without symmetrization
%           i.e. all possible values of (a,b,gamma,alpha,beta) 
%           s.t. (a+b \leq Nmax)
num_rows = size(S,1); 

% num_cols: Number of basis functions after symmetrization
%           expressed as a linear combination of all basis functions, 
%           i.e. all possible (a,b,gamma,alpha,beta) s.t. (a+b \leq Nmax)
num_cols = size(S,2);

% tot_inds: Indices corresponding to different (a,b,gamma,alpha,beta)
%           taken from the symm_orders array.
tot_inds = mbp_inds_ab_array(symm_orders);

if num_cols > 0
tot_Uprops = zeros(nnz(S),6);
ind_ranges = zeros(num_cols,2);

st1 = 1;
for ct1 = 1:num_cols
%     ct1
    ind1 = find(abs(S(:,ct1))>0);
    st2 = st1 - 1 + size(ind1,1);
    U_props = tot_inds(ind1,:);
    tot_Uprops(st1:st2,:) = U_props(:,[1,3:end]);
    ind_ranges(ct1,:) = [st1, st2];
    st1 = st2 + 1;
end
tot_Uprops(st2+1:end,:) = [];

mat_name = [data_fname0, ...
    'Sarr_MabInds_gbnull_zero_aPLUSb_max_',num2str(Nmax),'.mat'];
save(mat_name, 'tot_Uprops', 'ind_ranges','-v7.3');
end
end