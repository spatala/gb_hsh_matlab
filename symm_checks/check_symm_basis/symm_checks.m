function [] = symm_checks(pt_grp, Nmax, coeffs_typ, symm_typ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Steps:
%    1) Load a random grain boundary (g1, g2).
%    2) Create all the symmetrically equivalent grain boundaries
%    3) Compute MBP basis functions for all the symmetrically equivalent
%       grain boundaries.
%    4) Check that all the basis functions return the same vectors.
% 
% Input:
%    1) pt_grp
%    2) Nmax
%    3) coeffs_typ: string
%        'nmax' or 'aPLUSb_max'
%    4) symm_typ: string
%        'cryst', 'cryst_ges', 'zero', 'const'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
top_dir = get_top_dir();
util_dir = strcat(top_dir,'Util_functions','/');
addpath(genpath(util_dir));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%% Get data files
data_fname = [top_dir,'data_files/ptgrp_',pt_grp,'/'];

%%%% Load random GB parameter
s1 = load([top_dir,'data_files/GB_Parameters/rand_gb_rots.mat']); 
rot_mats = s1.rot_mats;
rots1 = rot_mats(:,:,floor(size(rot_mats,3)*rand()));
g1 = rots1(:,1:3); g2 = rots1(:,4:6);

%%%% Generate the symmetrically equivalent boundary parameters
%%%% The last parameter in get_symm_rots: opt
%%%%    1: Do not consider grain exchange symmetry
%%%%    2: Consider grain exchange symmetry
if strcmp(symm_typ, 'cryst')
    [symm_rots, ~] = get_symm_rots(g1,g2, pt_grp, data_fname,1);
else
    [symm_rots, ~] = get_symm_rots(g1,g2, pt_grp, data_fname,2);
end

s_angs = convert_gbrots(symm_rots);
SMvec=calc_Mvec_symm(top_dir, pt_grp, Nmax, ...
    coeffs_typ, symm_typ, s_angs);

disp(norm(full(SMvec - SMvec(1,:)))/size(SMvec,1))

rmpath(genpath(util_dir));

end