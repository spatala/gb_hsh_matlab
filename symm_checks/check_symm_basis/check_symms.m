%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check various symmetries of the basis functions
%
% This code will call the `symm_checks` and `symm_checks_nullgb_constraint` 
% for various symmetries, for a given `pt_grp`, `Nmax` and `coeffs_typ`
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;


pt_grp = 'Oh'; Nmax = 16; coeffs_typ = 'aPLUSb_max';

symm_typ = 'cryst';
symm_checks(pt_grp, Nmax, coeffs_typ, symm_typ)

symm_typ = 'cryst_ges';
symm_checks(pt_grp, Nmax, coeffs_typ, symm_typ)

symm_typ = 'zero';
symm_checks(pt_grp, Nmax, coeffs_typ, symm_typ)

symm_typ = 'const';
symm_checks(pt_grp, Nmax, coeffs_typ, symm_typ)

constraint = 'zero';
symm_checks_nullgb_constraint(pt_grp, Nmax, coeffs_typ, constraint)

constraint = 'const';
symm_checks_nullgb_constraint(pt_grp, Nmax, coeffs_typ, constraint)