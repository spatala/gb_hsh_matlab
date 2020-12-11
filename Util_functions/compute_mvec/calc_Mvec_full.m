function [Mvec, SMvec]=calc_Mvec_full(top_dir, pt_grp, Nmax, ...
    coeffs_typ, null_typ, rots)
% 
% Returns the full and symmetrized basis functions for the grain boundary 
% space using the SO3 X SO3 parameterization.
% 
% - Input:
%   + top_dir: Directory to find the symmetrized basis functions
%   + pt_grp: The crystallographic point group ('Oh')
%   + Nmax: all basis functions such that a+b \leq Nmax
%   + coeffs_typ: 'aPLUSb_max'
%   + null_typ: 'zero' or 'const'
%     - 'const': Crystal and GES and constant no-boundary condition
%     - 'null':  Crystal and GES and null no-boundary condition
%   + rots: gA and gB, orientations of the two crystals defining the GB.
%     - For N rotations, the array of size `N X 3 X 6`
% 
% - Output:
%   + Mvec: Full basis vectors
%   + SMvec: Symmetrized basis vectors
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];
data_fname0 = [data_fname,coeffs_typ,'_',num2str(Nmax),'/'];
mat_name = [data_fname0,'symm_ab_',pt_grp,'_',...
    coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name); symm_orders = s1.symm_orders;
mat_name = [data_fname0, 'Sarr_cryst_ges_gbnull_',...
    null_typ,'_',coeffs_typ,'_',num2str(Nmax),'.mat'];
s1 = load(mat_name); S = (s1.S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_val = symm_orders(:,1)'; 
b_val = symm_orders(:,2)'; 
c_val = min(a_val, b_val);
num_cols = sum((2*a_val+1).*(2*b_val+1).*(2*c_val+1));
tot_inds = mbp_inds_ab_array(symm_orders);

n_rots = size(rots,3);
ma2 = zeros(n_rots, 5);
for ct2=1:n_rots
    g2_1 = rots(:,1:3,ct2);
    g2_2 = rots(:,4:6,ct2);
    ma2(ct2,:) = rots_to_angs(g2_1, g2_2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsymm_evs = size(S,2);
nbfuncs = size(S,1);

SMvec = zeros(n_rots,nsymm_evs);
Mvec = zeros(n_rots,nbfuncs);

parfor ct2=1:n_rots
    Mvec2 = zeros(1,num_cols);
    for ct1 = 1:size(symm_orders,1)
        a=a_val(ct1);
        b = b_val(ct1);
        M2 = mbp_basis(a, b, [ma2(ct2,1), ma2(ct2,2), ...
            ma2(ct2,3), ma2(ct2,4), ma2(ct2,5)]);
        cond1 = (tot_inds(:,3)==a & tot_inds(:,4)==b);
        tinds1 = find(cond1); 
        ind_start = tinds1(1); 
        ind_stop = tinds1(end);
        Mvec2(ind_start:ind_stop) = M2;
    end
    Mvec(ct2,:) = Mvec2;
    SMvec(ct2,:) = Mvec2*S;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
