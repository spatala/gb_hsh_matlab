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
1. `clean.m` : Attempts to minimize the number of nonzero coefficients of the eigenvectors in `S` while maintaining orthonormality.
2. `clebsch_gordan.m`
3. `spnull.m`

## param_conversions
1. `check_conv.m`
2. `mbp_to_rots.m`: Convert GB parameters from Misorientation-Boundary plane (MBP) to SO3xSO3 parametrization (rots).
3. `rotmat_to_angs.m`
4. `rots_to_angs.m`
5. `rots_to_mbp.m`


## rotation_Umat

This folder contains functions to construct SO(3) irreducible representatives (matrices).

### Functions
1. `angles_to_CK.m`:returns `U^{1/2}` for the rotation given by a rotation angle and axis. Rows and cols ordered by increasing `m'` and `m`.
2. `CK_to_angles.m`: returns the rotation angle and axis for the given `U^{1/2}` with rows and cols ordered by increasing m' and m.
3. `rotation_mat.m`: returns the rotation matrix `U` of the order `j`, with a rotation angle `w` and rotation axis given by polar and azimuthal angles `th` and `ph`.
4. `wigner_little_d.m`: returns the Wigner `d` matrix for the given order and rotation angle, with rows and columns labeled in increasing values of `m`. The second argument is an estimate of the magnitude of numerical error.

### Functions To be modified or removed:
1. **To be modified**
	+ `rotation_wo_svd.m`: Function to compute elements of SO(3) Irreducible representative  `U^{a_val}_{alp_val, al_val}` for a rotation given by (q,Q) quaternion.
2. **To be removed**
	+ `rotation.m`: Constructs irreducible representation of SO(3) of dimension (N+1) for the rotation specified by axis `ax` and angle `an`.





# Functions List

## Util_functions

### To Document

1. `spnull`
2. `sp_null`
3. `rotation`

### Doucmented

`get_symmgen_mats`
`get_symm_rots`
`check_conv`
`mbp_inds_ab_array`
`mbp_ab`
`mbp_basis`
`convert_gbrots`
`compute_Mvec`
`clean`
`clebsch_gordan`
`sp_orth`
`mbp_to_rots`
`rotmat_to_angs`
`rots_to_angs`
`rots_to_mbp`
`angles_to_CK`
`CK_to_angles`
`rotation_mat`
`rotation_wo_svd`
`wigner_little_d`
`get_symmgen_angs`

 












