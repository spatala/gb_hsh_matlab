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
1. `CK_to_angles.m`
2. `angles_to_CK.m`
3. `rotation.m`: Constructs irreducible representation of SO(3) of dimension (N+1) for the rotation specified by axis ax and angle an.
4. `rotation_mat.m`
5. `rotation_wo_svd.m`:Function to compute elements of SO(3) Irreducible representative  `U^{a_val}_{alp_val, al_val}` for a rotation given by (q,Q) quaternion.
6. `wigner_little_d.m`















