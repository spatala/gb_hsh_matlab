# Code Descriptions

## `calc_Mfunc_Rab(Rab, a,b)`

Get the value of MBP functions for a given `(a,b)` and `R^(ab)` matrix.

## `mbp_inds(N)`

Get indices for MBP functions all the way to a=b=N.

## `get_symmgen_mats(pt_grp)`

Function to set the generators for each point group (only proper rotation elements are set.)

## `mbp_ab(a,b)`

Get indices for MBP functions with fixed a and b.

## `mbp_inds_ab_array(symm_orders)`

Get indices for MBP function with values of (a,b) given by `symm_orders`.

## `mbp_to_rots(mbp)`

Convert GB parameters from Misorientation-Boundary plane (MBP) to SO3xSO3 parametrization (rots).

## `so4_irrep(g1,g2,Na,Nb)`

Function to compute SO(4) Irreducible representative with g1, g2 rotation matrices. The output is a square matrix with (Na+1)(Nb+1) rows.

## `rotation_wo_svd(a,lb,q, Q,a_val, alp_val, al_val)`

Function to compute elements of SO(3) Irreducible representative  `U^{a_val}_{alp_val, al_val}` for a rotation given by (q,Q) quaternion.

## `rotation(ax,an,N)`

Constructs irreducible representation of SO(3) of dimension (N+1) for the rotation specified by axis ax and angle an. 





