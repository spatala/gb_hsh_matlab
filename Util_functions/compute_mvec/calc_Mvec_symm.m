function SMvec=calc_Mvec_symm(top_dir, pt_grp, Nmax, ...
    coeffs_typ, symm_typ, s_angs)
% 
% Computes symmetrized MBP basis functions for given set of
% (a, b, gamma, alpha, beta) for a given order (a+b \leq Nmax)
% of symmetrized MBP basis function.
% 
% - Input:
%   + top_dir: Directory to find the symmetrized basis functions
%   + pt_grp: The crystallographic point group ('Oh')
%   + Nmax: all basis functions such that a+b \leq Nmax
%   + coeffs_typ: 'aPLUSb_max'
%   + symm_typ: string ('cryst' or 'cryst_ges' or 'const' or 'zero')
%     - specify which level of symmetry to use to compute basis functions
%     - 'cryst': Only crystal point group symmetries
%     - 'cryst_ges': Crystal and GES
%     - 'const': Crystal and GES and constant no-boundary condition
%     - 'null':  Crystal and GES and null no-boundary condition
%   + s_angs: struct of N(=nrots)-arrays (8 elements)
%     - ((s.a1, s.lb1, s.q1, s.Q1), (s.a2, s.lb2, s.q2, s.Q2))
% 
% - Output:
%   + SMvec: symmetrized basis functions
%     - The rows coresspond to different rotations.
%     - The columns correspond to the number of symmetrized basis functions.
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,coeffs_typ,'_',num2str(Nmax),'/'];
mat_name = [data_fname0,'symm_ab_',pt_grp,'_',...
    coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;

if strcmp(symm_typ,'cryst')
    mat_name = [data_fname0, 'Sarr_abc_combined_csymm_',...
        coeffs_typ,'_',num2str(Nmax),'.mat'];
    mat_name1 = [top_dir, 'data_files/','ptgrp_',pt_grp,'/', coeffs_typ,'_',...
    num2str(Nmax),'/','Sarr_MabInds_',symm_typ,'_',...
    coeffs_typ,'_',num2str(Nmax),'.mat'];
elseif strcmp(symm_typ,'cryst_ges')
    mat_name = [data_fname0, 'Sarr_cryst_ges_',...
        coeffs_typ,'_',num2str(Nmax),'.mat'];
    mat_name1 = [top_dir, 'data_files/','ptgrp_',pt_grp,'/', coeffs_typ,'_',...
    num2str(Nmax),'/','Sarr_MabInds_',symm_typ,'_',...
    coeffs_typ,'_',num2str(Nmax),'.mat'];
else
    mat_name = [data_fname0, 'Sarr_cryst_ges_gbnull_',symm_typ,'_',...
        coeffs_typ,'_',num2str(Nmax),'.mat'];
    mat_name1 = [top_dir, 'data_files/','ptgrp_',pt_grp,'/', coeffs_typ,'_',...
    num2str(Nmax),'/','Sarr_MabInds_gbnull_',symm_typ,'_',...
    coeffs_typ,'_',num2str(Nmax),'.mat'];
end

s1 = load(mat_name); S = (s1.S);
s1 = load(mat_name1); tot_Uprops = s1.tot_Uprops;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tot_inds = mbp_inds_ab_array(symm_orders);
vec_inds = tot_inds(unique(tot_Uprops(:,1)),:);
nrots = size(s_angs.a1,1);
% s = convert_gbrots(rots);
num_rows = size(S,1);
Mvec = compute_Mvec(s_angs, nrots, num_rows, vec_inds);
SMvec = Mvec*S;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
