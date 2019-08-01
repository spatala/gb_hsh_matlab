# Code Descriptions

## compute_mvec
1. `mbp_basis.m`
2. `convert_gbrots.m`
3. `compute_Mvec.m`
4. `compute_symm_Mvec.m`


## symmetry_generators
1. `get_symm_rots.m`
2. `get_symmgen_angs.m`: Function to set the generators (as angles - (w, th, ph)) for each point group (only proper rotation elements are set.)
3. `get_symmgen_mats.m`: Function to set the generators (as matrices) for each point group (only proper rotation elements are set.)

## mbp_indices
1. `mbp_ab.m`: Get indices for MBP functions with fixed a and b.
2. `mbp_inds_array.m`: Get indices for MBP function with values of (a,b) given by `symm_orders`.

## misc_funcs
1. `clean.m` : Attempts to minimize the number of nonzero coefficients of the eigenvectors in S while maintaining orthonormality.
2. `clebsch_gordan.m`
3. `spnull.m`

## param_conversions
1. `check_conv.m`
2. `mbp_to_rots.m`: Convert GB parameters from Misorientation-Boundary plane (MBP) to SO3xSO3 parametrization (rots).
3. `rotmat_to_angs.m`
4. `rots_to_angs.m`
5. `rots_to_mbp.m`


## rotation_Umat
+ `angles_to_CK.m`:returns `U^{1/2}` for the rotation given by a rotation angle and axis. Rows and cols ordered by increasing `m'` and `m`.
+ `CK_to_angles.m`: returns the rotation angle and axis for the given `U^{1/2}` with rows and cols ordered by increasing m' and m.
+ `rotation_mat.m`: returns the rotation matrix `U` of the order `j`, with a rotation angle `w` and rotation axis given by polar and azimuthal angles `th` and `ph`.
+ `wigner_little_d.m`: returns the Wigner `d` matrix for the given order and rotation angle, with rows and columns labeled in increasing values of `m`. The second argument is an estimate of the magnitude of numerical error.
### To be modified or deprecated:
+ `rotation.m`: Constructs irreducible representation of SO(3) of dimension (N+1) for the rotation specified by axis ax and angle an. **(Should be deleted in future versions.)**
+ `rotation_wo_svd.m`:Function to compute elements of SO(3) Irreducible representative  `U^{a_val}_{alp_val, al_val}` for a rotation given by (q,Q) quaternion. **(Has to be modified!)**

















