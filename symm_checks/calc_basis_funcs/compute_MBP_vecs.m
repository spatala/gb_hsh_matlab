%%%%%%%%%
% Codes to illustrate the computation of MBP basis vectors
% 1) Compute full basis vectors and symmetrized using mbp_basis
% 2) Compute symmetrized basis vectors using compute_Mvec
%
clear all; clc;

top_dir = get_top_dir();
du1 = [top_dir, 'Util_functions/'];
addpath(genpath(du1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load random GB parameters
s1 = load([top_dir,'data_files/GB_Parameters/rand_gb_rots.mat']); 
rots = s1.rot_mats;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pt_grp = 'Oh'; Nmax = 16; coeffs_typ = 'aPLUSb_max'; null_typ = 'zero';

tic;
[Mvec, SMvec]=calc_Mvec_full(top_dir, pt_grp, Nmax, coeffs_typ, null_typ, rots);
toc;

tic;
s_angs = convert_gbrots(rots);
SMvec1=calc_Mvec_symm(top_dir, pt_grp, Nmax, coeffs_typ, null_typ, s_angs);
toc;

display(max(max(abs(SMvec - SMvec1))));
rmpath(genpath(du1));